% bpsk_hf.m
% David Rowe Mar 2017
%
% Rate Rs BPSK simulation to explore phase over multiple carriers in HF channel

1;

function sim_out = run_sim(sim_in)
  Rs = 100;

  Nbits = sim_in.Nbits;
  EbNodB = sim_in.EbNodB;
  verbose = sim_in.verbose;
  hf_en = sim_in.hf_en;
  Nc = sim_in.Nc;          % Number of cols, aka number of carriers
  Nr = Nbits/Nc;           % Number of rows to get Nbits total

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

    rx = tx;

    if hf_en

      % simplified rate Rs simulation model that doesn't include
      % ISI, just freq filtering.  We assume perfect phase estimation
      % so it's just amplitude distortion.
      
      % Note Rs carrier spacing, sample rate is Rs

      hf_model = zeros(Nr,Nc); phase_est = zeros(1,Nr);
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
          for c=1:Nc
            rx(i,c) = rx(i,c) .* hf_model(i,c) * exp(-j*phase_est(i,c));
          end

        else

          % just apply amplitude fading, assume we have ideal phase est

          for c=1:Nc
            rx(i,c) = rx(i,c) .* abs(hf_model(i,c));
          end
        end
      end
    end

    rx += noise;
    rx_bits = real(rx) < 0;
    errors = xor(tx_bits, rx_bits);
    Nerrs = sum(errors);
    printf("EbNodB: %3.2f BER: %4.3f Nbits: %d Nerrs: %d\n", EbNodB(nn), sum(Nerrs)/Nbits, Nbits, sum(Nerrs));

    if verbose
      figure(1); clf; 
      plot(rx,'+');
      axis([-2 2 -2 2]);
      figure(2); clf; 
      plot(abs(hf_model));
      figure(3); clf; 
      plot(angle(hf_model));
      hold on; plot(phase_est,'+'); hold off;
    end

    sim_out.ber(nn) = sum(Nerrs)/Nbits;
  end
endfunction


function run_curves
  sim_in.verbose = 0;
  sim_in.Nbits = 12000;
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
  sim_in.Nbits = 3000;
  sim_in.EbNodB = 6;
  sim_in.verbose = 0;
  sim_in.hf_en = 1;
  sim_in.Nc = 3;
  sim_in.av_phase = 1;

  run_sim(sim_in);
end


format;
more off;

run_single
%run_curves



