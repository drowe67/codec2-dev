% bpsk_hf.m
% David Rowe Mar 2017
%
% Rate Rs BPSK simulation to explore phase estimation
% over multiple carriers in HF channel

#{
 TODO:
   [X] sim pilot based phase est using known symbols
   [X] test AWGN BER with averaging pilots from adj carriers
   [X] refactor to insert pilot rows
   [X] add border cols, not used for data
   [X] centre est on current carrier, extend to > 3
   [X] test single points
       + 1dB IL @ 6dB HF, 0.4 dB @ 2dB AWGN
   [ ] try linear interpolation
   [ ] try longer time windows
   [ ] try combining mod stripping phase est inside frame
   [ ] curves taking into account pilot losses
   [ ] remove border carriers, interpolate edge carrier
#}

1;

function sim_out = run_sim(sim_in)
  Rs = 100;

  Nbits = sim_in.Nbits;
  EbNodB = sim_in.EbNodB;
  verbose = sim_in.verbose;
  hf_en = sim_in.hf_en;
  hf_phase = sim_in.hf_phase;
  Ns = sim_in.Ns;          % step size for pilots
  Nc = sim_in.Nc;          % Number of cols, aka number of carriers
  Nr = Nbits/Nc;           % Number of rows to get Nbits total
  phase_offset = sim_in.phase_offset;

  if verbose
    printf("Nbits: %d\n", Nbits);
    printf("Nc: %d\n", Nc);
    printf("Nr: %d\n", Nr);
    printf("Ns: %d (step size for pilots, Ns-1 data symbols between pilots)\n", Ns);
  end

  % check if Nbits fit neatly into carriers

  assert(Nbits/Nc == floor(Nbits/Nc), "Nbits/Nc must be an integer");

  % check if bits fit neatly into frames with rows of pilots above and below
  % PPP
  % DDD
  % DDD
  % PPP
  
  Nbitsperframe = (Ns-1)*Nc;
  printf("Nbitsperframe: %d\n", Nbitsperframe);
  Nframes = Nbits/Nbitsperframe;
  Nrowsperframe = Nbitsperframe/Nc;
  printf("Nrowsperframe: %d\n", Nrowsperframe);

  % check if Nbits fit neatly into frames delineated by pilots

  assert(Nframes == floor(Nframes), "Nbits/Nbits/frame must be an integer");
  printf("Nframes: %d\n", Nframes);

  Nrp = Nr + Nframes + 1;  % number of rows once pilots inserted
                           % extra row of pilots at end
  printf("Nrp: %d (number of rows including pilots)\n", Nrp);

  % set up HF model

  if hf_en

    % some typical values

    dopplerSpreadHz = 1.0; path_delay = 1E-3*Rs;

    spread1 = doppler_spread(dopplerSpreadHz, Rs, Nrp+10);
    spread2 = doppler_spread(dopplerSpreadHz, Rs, Nrp+10);

    % sometimes doppler_spread() doesn't return exactly the number of samples we need
 
    assert(length(spread1) >= Nrp, "not enough doppler spreading samples");
    assert(length(spread2) >= Nrp, "not enough doppler spreading samples");

    hf_gain = 1.0/sqrt(var(spread1)+var(spread2));
    % printf("nsymb: %d lspread1: %d\n", nsymb, length(spread1));
  end
  
  % construct an artificial phase countour for testing, linear across freq and time

  if sim_in.phase_test
    phase_test = ones(Nrp, Nc+2);
    for r=1:Nrp
      for c=1:Nc+2
        phase_test(r,c) = -pi/2 + c*pi/(Nc+2) + r*0.01*2*pi;
        phase_test(r,c) = phase_test(r,c) - 2*pi*floor((phase_test(r,c)+pi)/(2*pi));
      end
    end
  end

  % simulate for each Eb/No point ------------------------------------

  for nn=1:length(EbNodB)
    rand('seed',1);
    randn('seed',1);

    EbNo = 10 .^ (EbNodB(nn)/10);
    variance = 1/(EbNo/2);
    noise = sqrt(variance)*(0.5*randn(Nrp,Nc+2) + j*0.5*randn(Nrp,Nc+2));

    % generate tx bits and insert pilot rows and border cols

    tx_bits = []; tx_bits_np = [];
    for f=1:Nframes
      tx_bits = [tx_bits; ones(1,Nc+2)];
      for r=1:Nrowsperframe
        arowofbits = rand(1,Nc) > 0.5;
        tx_bits = [tx_bits; [1 arowofbits 1]];
        tx_bits_np = [tx_bits_np; arowofbits];
      end
    end      
    tx_bits = [tx_bits; [1 ones(1,Nc) 1]]; % final row of pilots
    [nr nc] = size(tx_bits);
    assert(nr == Nrp);
    %tx_bits

    tx = 2*tx_bits - 1;

    rx = tx * exp(j*phase_offset);

    if sim_in.phase_test
      rx = rx .* exp(j*phase_test);
    end

    if hf_en

      % simplified rate Rs simulation model that doesn't include
      % ISI, just freq filtering.
      
      % Note Rs carrier spacing, sample rate is Rs

      hf_model = zeros(Nr,Nc+2); phase_est = zeros(Nr,Nc);
      for r=1:Nrp
        for c=1:Nc+2
          w = 2*pi*c*Rs/Rs; 
          hf_model(r,c) = hf_gain*(spread1(r) + exp(-j*w*path_delay)*spread2(r));
        end
        
        if hf_phase
          rx(r,:) = rx(r,:) .* hf_model(r,:);
        else
          rx(r,:) = rx(r,:) .* abs(hf_model(r,:));
        end
      end
    end

    rx += noise;

    % pilot based phase est, we use known tx symbols as pilots ----------

    rx_corr = rx;

    if sim_in.pilot_phase_est

      % est phase from pilots either side of data symbols
      % adjust phase of data symbol
      % demodulate and count errors of just data
 
      phase_est_pilot_log = 10*ones(Nrp,Nc+2);
      phase_est_stripped_log = 10*ones(Nrp,Nc+2);
      phase_est_log = 10*ones(Nrp,Nc+2);
      for c=2:Nc+1
        for r=1:Ns:Nrp-Ns

          % estimate phase using average of 6 pilots in a rect 2D window centred
          % on this carrier
          % PPP
          % DDD
          % DDD
          % PPP
          
          cr = c-1:c+1;
          aphase_est_pilot_rect = sum(rx(r,cr)*tx(r,cr)') +  sum(rx(r+Ns,cr)*tx(r+Ns,cr)');

          % optionally use next step of pilots in past and future

          if sim_in.pilot_wide
            if r > Ns+1
              aphase_est_pilot_rect += sum(rx(r-Ns,cr)*tx(r-Ns,cr)');
            end
            if r < Nrp - 2*Ns
              aphase_est_pilot_rect += sum(rx(r+2*Ns,cr)*tx(r+2*Ns,cr)');
            end
          end

          aphase_est_pilot = angle(aphase_est_pilot_rect);

          % correct phase offset using phase estimate

          for rr=r+1:r+Ns-1
            phase_est_pilot_log(rr,c) = aphase_est_pilot;
            rx_corr(rr,c) = rx(rr,c) * exp(-j*aphase_est_pilot);
          end

          if sim_in.stripped_phase_est
            % Optional modulation stripping feed fwd phase estimation, to refine
            % pilot-based phase estimate.  Doing it after pilot based phase estimation
            % means we don't need to deal with ambiguity, which is difficult to handle
            % in low SNR channels.

            % Use vector of 7 symbols around current data symbol.  We could use a 2D
            % window if we can work out how best to correct with pilot-est and avoid
            % ambiguities

            for rr=r+1:r+Ns-1

              % extract a matrix of nearby samples with pilot-based offset removed
 
              amatrix = rx(max(1,rr-3):min(Nrp,rr+3),c) .* exp(-j*aphase_est_pilot);

              % modulation strip and est phase

              stripped = abs(amatrix) .* exp(j*2*angle(amatrix));
              aphase_est_stripped = angle(sum(sum(stripped)))/2;
              phase_est_stripped_log(rr,c) = aphase_est_stripped;

              % correct rx symbols based on both phase ests

              phase_est_log(rr,c) = angle(exp(j*(aphase_est_pilot+aphase_est_stripped)));
              rx_corr(rr,c) = rx(rr,c) * exp(-j*phase_est_log(rr,c));
            end       
          end % sim_in.stripped_phase_est

        end % r=1:Ns:Nrp-Ns

      end % c=2:Nc+1
    end % sim_in.pilot_phase_est

    % remove pilots to give us just data symbols
      
    rx_np = [];
    for r=1:Nrp
      if mod(r-1,Ns) != 0
        rx_np = [rx_np; rx_corr(r,2:Nc+1)];
      end
    end

    %phase_test
    %phase_est_log

    % calculate BER stats as a block, after pilots extracted

    rx_bits_np = real(rx_np) > 0;
    errors = xor(tx_bits_np, rx_bits_np);
    Nerrs = sum(sum(errors));

    printf("EbNodB: %3.2f BER: %4.3f Nbits: %d Nerrs: %d\n", EbNodB(nn), Nerrs/Nbits, Nbits, Nerrs);

    if verbose
      figure(1); clf; 
      plot(rx_np,'+');
      axis([-2 2 -2 2]);

      if hf_en
        figure(2); clf; 
        plot(abs(hf_model(:,2:Nc+1)));
      end

      if sim_in.pilot_phase_est
        figure(3); clf;
        plot(phase_est_log(:,2:Nc+1),'+', 'markersize', 10); 
        hold on; 
        plot(phase_est_pilot_log(:,2:Nc+1),'g+', 'markersize', 5); 
        if sim_in.stripped_phase_est
          plot(phase_est_stripped_log(:,2:Nc+1),'ro', 'markersize', 5); 
        end
        if sim_in.hf_phase
          plot(angle(hf_model(:,2:Nc+1)));
        end
        if sim_in.phase_test
          plot(phase_test(:,2:Nc+1));
        end
        axis([1 Nrp -pi pi]);
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
  sim_in.Nc = 7;
  sim_in.Ns = 8;
  sim_in.Nbits = 1000*sim_in.Nc*(sim_in.Ns-1);
  sim_in.EbNodB = 6;
  sim_in.verbose = 1;
  sim_in.pilot_phase_est = 0;
  sim_in.pilot_wide = 1;
  sim_in.stripped_phase_est = 0;
  sim_in.phase_offset = 0;
  sim_in.phase_test = 0;
  sim_in.hf_en = 1;
  sim_in.hf_phase = 0;

  run_sim(sim_in);
end


format;
more off;

run_single
%run_curves



