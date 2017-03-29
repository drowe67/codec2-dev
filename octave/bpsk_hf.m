% bpsk_hf.m
% David Rowe Mar 2017
%
% Rate Rs BPSK simulation to explore phase estimation
% over multiple carriers in HF channel

#{
 TODO:
   [ ] sim pilot based phase est using known symbols
   [ ] test AWGN BER with averaging pilots from adj carriers
#}

1;

function sim_out = run_sim(sim_in)
  Rs = 100;

  Nbits = sim_in.Nbits;
  EbNodB = sim_in.EbNodB;
  verbose = sim_in.verbose;
  hf_en = sim_in.hf_en;
  Ns = sim_in.Ns;          % step size for pilots
  Nc = sim_in.Nc;          % Number of cols, aka number of carriers
  Nr = Nbits/Nc;           % Number of rows to get Nbits total
  phase_offset = sim_in.phase_offset;

  assert(Nbits/Nc == floor(Nbits/Nc), "Nbits/Nc must be an integer");
  
  % set up HF model

  if hf_en

    % some typical values

    dopplerSpreadHz = 1.0; path_delay = 1E-3*Rs;

    spread1 = doppler_spread(dopplerSpreadHz, Rs, Nr*2);
    spread2 = doppler_spread(dopplerSpreadHz, Rs, Nr*2);
    hf_gain = 1.0/sqrt(var(spread1)+var(spread2));
    % printf("nsymb: %d lspread1: %d\n", nsymb, length(spread1));
  end
  
  % simulate for each Eb/No point

  for nn=1:length(EbNodB)
    rand('seed',1);
    randn('seed',1);

    EbNo = 10 .^ (EbNodB(nn)/10);
    variance = 1/(EbNo/2);
    noise = sqrt(variance)*(0.5*randn(Nr,Nc) + j*0.5*randn(Nr,Nc));

    tx_bits = rand(Nr,Nc) > 0.5;
    tx = 1 - 2*tx_bits;

    rx = tx * exp(j*phase_offset);

    if hf_en

      % simplified rate Rs simulation model that doesn't include
      % ISI, just freq filtering.  We assume perfect phase estimation
      % so it's just amplitude distortion.
      
      % Note Rs carrier spacing, sample rate is Rs

      hf_model = zeros(Nr,Nc); phase_est = zeros(Nr,Nc);
      for i=1:Nr
        for c=1:Nc
          w = 2*pi*c*Rs/Rs; 
          hf_model(i,c) = hf_gain*(spread1(i) + exp(-j*w*path_delay)*spread2(i));
        end

        if sim_in.av_phase

          % Simulate obtaining phase est by averaging hf channel model
          % samples at carrier frequencies adjacent to each carrier.
          % In case of edge carriers we sample HF channel model to one
          % side of carrier for now, for a real world estimate we'd
          % need to interpolate phase to edge carriers somehow

          for c=1:Nc
            phase_sum = 0;
            for cc = c-1:c+1
              w = 2*pi*cc*Rs/Rs; 
              ahf_model = hf_gain*(spread1(i) + exp(-j*w*path_delay)*spread2(i));
              phase_sum += ahf_model;
            end
            phase_est(i,c) = angle(phase_sum); 
          end
          rx(i,:) = rx(i,:) .* hf_model(i,:) .* exp(-j*phase_est(i,:));

        else

          % just apply amplitude fading, assume we have ideal phase est

          rx(i,:) = rx(i,:) .* abs(hf_model(i,:));
        end
      end
    end

    rx += noise;

    Nbits_np = Nbits;
    rx_np = rx;
    tx_bits_np = tx_bits;

    % pilot based phase est, we use known tx symbols as pilots ----------

    if sim_in.pilot_phase_est

      % est phase from pilots either side of data symbols
      % adjust phase of data symbol
      % demodulate and count errors of just data
 
      phase_est_log = zeros(Nr,Nc);
      for c=1:Nc
        for r=1:Ns:Nr-Ns

          % estimate phase using average of two pilots at either end of frame

          %aphase_est_rect = rx(r,c)*tx(r,c)' +  rx(r+Ns,c)*tx(r+Ns,c)';
          aphase_est_rect = sum(rx(r,:)*tx(r,:)') +  sum(rx(r+Ns,:)*tx(r+Ns,:)');
          aphase_est = angle(aphase_est_rect);

          % apply this phase estimate to symbols in frame, and check for error

          for rr=r+1:r+Ns-1
            rx(rr,c) *= exp(-j*aphase_est);
            phase_est_log(rr,c) = aphase_est;
          end
        end
      end

      % build up tx_bits and rx symbols matrices without pilots
      
      tx_bits_np = []; rx_np = []; Nbits_np = 0;
      for r=1:Nr
        if mod(r-1,Ns) != 0
          tx_bits_np = [tx_bits_np; tx_bits(r,:)];
          rx_np = [rx_np; rx(r,:)];
          Nbits_np += Nc;
        end
      end
    end

    % calculate BER stats as a block, after pilots extracted

    rx_bits_np = real(rx_np) < 0;
    errors = xor(tx_bits_np, rx_bits_np);
    Nerrs = sum(sum(errors));

    printf("EbNodB: %3.2f BER: %4.3f Nbits: %d Nerrs: %d\n", EbNodB(nn), Nerrs/Nbits_np, Nbits_np, Nerrs);

    if verbose
      figure(1); clf; 
      plot(rx_np,'+');
      axis([-2 2 -2 2]);
      if sim_in.pilot_phase_est
        figure(2); clf;
        plot(phase_est_log,'+');
      end
      if hf_en
        figure(2); clf; 
        plot(abs(hf_model));
        figure(3); clf; 
        plot(angle(hf_model));
        hold on; plot(phase_est,'+'); hold off;
      end
    end

    sim_out.ber(nn) = sum(Nerrs)/Nbits;
  end
endfunction


function run_curves
  sim_in.verbose = 0;
  sim_in.Nbits = 90000;
  sim_in.EbNodB = 2:8;
  sim_in.hf_en = 1;
  sim_in.Nc = 3;

  sim_in.av_phase = 0;
  bpsk_hf = run_sim(sim_in);

  sim_in.av_phase = 1;
  bpsk_hf_av_phase = run_sim(sim_in);

  figure(3); clf;
  semilogy(sim_in.EbNodB, bpsk_hf.ber,'b+-;BPSK HF;');
  hold on;
  semilogy(sim_in.EbNodB, bpsk_hf_av_phase.ber,'g+-;BPSK HF average phase;');
  hold off;
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid;
  legend('boxoff');
end


function run_single
  data_rows = 1000;
  sim_in.Nc = 3;
  sim_in.Ns = 6;
  sim_in.Nbits = (data_rows*sim_in.Ns+1)*sim_in.Nc;
  sim_in.EbNodB = 2;
  sim_in.verbose = 1;
  sim_in.pilot_phase_est = 1;
  sim_in.phase_offset = -pi;
  sim_in.hf_en = 0;
  sim_in.av_phase = 0;

  run_sim(sim_in);
end


format;
more off;

run_single
%run_curves



