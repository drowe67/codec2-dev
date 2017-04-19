% bpsk_hf_fs.m
% David Rowe Mar 2017
%
% Rate Fs BPSK simulation, development of bpsk_hf_rs.m

#{
 TODO:
   [X] strip back experimental stuff to just features we need
   [X] ZOH/integrator
   [X] OFDM up and down conversion
   [X] rate Fs HF model and HF results
   [X] add QPSK
   [X] add CP
   [ ] adjust waveform parameters for real world
   [ ] Nsec run time take into account CP
   [ ] plot phase ests
   [ ] handle border carriers
       [ ] start with phantom carriers 
           + but unhappy with 1800Hz bandwidth
       [ ] also try interpolation or just single row
   [ ] compute SNR and PAPR
   [ ] timing estimator
   [ ] acquisition & freq offset estimation
   [ ] SSB bandpass filtering
#}

1;

% Gray coded QPSK modulation function

function symbol = qpsk_mod(two_bits)
    two_bits_decimal = sum(two_bits .* [2 1]); 
    switch(two_bits_decimal)
        case (0) symbol =  1;
        case (1) symbol =  j;
        case (2) symbol = -j;
        case (3) symbol = -1;
    endswitch
endfunction


% Gray coded QPSK demodulation function

function two_bits = qpsk_demod(symbol)
    bit0 = real(symbol*exp(j*pi/4)) < 0;
    bit1 = imag(symbol*exp(j*pi/4)) < 0;
    two_bits = [bit1 bit0];
endfunction


% Returns the most likely place for the start of a frame, as a
% a candidate for coarse frame sync.  Combines two frames pilots
% so we need at least Nsamperframe+M+Ncp samples in rx

function [ct_est foff_est] = coarse_sync(states, rx, rate_fs_pilot_samples)
    Nsamperframe = states.Nsamperframe; M = states.M; Ncp = states.Ncp; Fs = states.Fs;
    verbose  = states.verbose;

    Ncorr = length(rx) - (Nsamperframe+M+Ncp) + 1;
    assert(Ncorr > 0);
    corr = zeros(1,Ncorr);
    for i=1:Ncorr
      corr(i)   = abs(rx(i:i+M+Ncp-1) * rate_fs_pilot_samples');
      corr(i)  += abs(rx(i+Nsamperframe:i+Nsamperframe+M+Ncp-1) * rate_fs_pilot_samples');
    end

    [mx ct_est] = max(abs(corr));

    C  = abs(fft(rx(ct_est:ct_est+M+Ncp-1) .* conj(rate_fs_pilot_samples), Fs));
    C += abs(fft(rx(ct_est+Nsamperframe:ct_est+Nsamperframe+M+Ncp-1) .* conj(rate_fs_pilot_samples), Fs));

    fmax = 30;
    [mx_pos foff_est_pos] = max(C(1:fmax));
    [mx_neg foff_est_neg] = max(C(Fs-fmax+1:Fs));

    if mx_pos > mx_neg
      foff_est = foff_est_pos - 1;
    else
      foff_est = foff_est_neg - fmax - 1;
    end

    if verbose
      %printf("ct_est: %d\n", ct_est);
      figure(6); clf;
      plot(abs(corr))
      figure(7)
      plot(C)
      axis([0 200 0 0.4])
    end
endfunction


function [sim_out rate_fs_pilot_samples rx] = run_sim(sim_in)
  Fs = 8000;
  Rs = sim_in.Rs;
  Tcp = sim_in.Tcp;
  M  = Fs/Rs;
  Ncp = Tcp*Fs;
  woffset = 2*pi*sim_in.foff_hz/Fs;

  bps = sim_in.bps;
  EbNodB  = sim_in.EbNodB;
  verbose = sim_in.verbose;
  hf_en   = sim_in.hf_en;
  timing_en = sim_in.timing_en;

  Ns = sim_in.Ns;          % step size for pilots
  Nc = sim_in.Nc;          % Number of cols, aka number of carriers

  Nbitsperframe = (Ns-1)*Nc*bps;
  Nrowsperframe = Nbitsperframe/(Nc*bps);
  if verbose
    printf("Rs:.........: %4.2f\n", Rs);
    printf("M:..........: %d\n", M);
    printf("bps:.........: %d\n", bps);
    printf("Nbitsperframe: %d\n", Nbitsperframe);
    printf("Nrowsperframe: %d\n", Nrowsperframe);
  end

  % Important to define run time in seconds so HF model will evolve the same way
  % for different pilot insertion rates.  So lets work backwards from approx
  % seconds in run to get Nbits, the total number of payload data bits

  % frame has Ns-1 data symbols between pilots, e.g. for Ns=3: 
  %
  % PPP
  % DDD
  % DDD
  % PPP

  Nrows = sim_in.Nsec*Rs;
  Nframes = floor((Nrows-1)/Ns);
  Nbits = Nframes * Nbitsperframe;    % number of payload data bits

  Nr = Nbits/(Nc*bps);                % Number of data rows to get Nbits total

  % double check if Nbits fit neatly into carriers

  assert(Nbits/(Nc*bps) == floor(Nbits/(Nc*bps)), "Nbits/(Nc*bps) must be an integer");

  Nrp = Nr + Nframes + 1;  % number of rows once pilots inserted
                           % extra row of pilots at end
  Nsamperframe =  (Nrowsperframe+1)*(M+Ncp);

  if verbose
    printf("Nc.....: %d\n", Nc);
    printf("Ns.....: %d (step size for pilots, Ns-1 data symbols between pilots)\n", Ns);
    printf("Nr.....: %d\n", Nr);
    printf("Nbits..: %d\n", Nbits);
    printf("Nframes: %d\n", Nframes);
    printf("Nrp....: %d (number of rows including pilots)\n", Nrp);
  end

  % generate same pilots each time

  rand('seed',1);
  pilots = 1 - 2*(rand(1,Nc+2) > 0.5);

  % set up HF model ---------------------------------------------------------------

  if hf_en

    % some typical values, or replace with user supplied

    dopplerSpreadHz = 1.0; path_delay_ms = 1;

    if isfield(sim_in, "dopplerSpreadHz") 
      dopplerSpreadHz = sim_in.dopplerSpreadHz;
    end
    if isfield(sim_in, "path_delay_ms") 
      path_delay_ms = sim_in.path_delay_ms;
    end
    path_delay_samples = path_delay_ms*Fs/1000;
    printf("Doppler Spread: %3.2f Hz Path Delay: %3.2f ms %d samples\n", dopplerSpreadHz, path_delay_ms, path_delay_samples);

    % generate same fading pattern for every run

    randn('seed',1);

    spread1 = doppler_spread(dopplerSpreadHz, Fs, (sim_in.Nsec*(M+Ncp)/M+0.2)*Fs);
    spread2 = doppler_spread(dopplerSpreadHz, Fs, (sim_in.Nsec*(M+Ncp)/M+0.2)*Fs);

    % sometimes doppler_spread() doesn't return exactly the number of samples we need
 
    assert(length(spread1) >= Nrp*M, "not enough doppler spreading samples");
    assert(length(spread2) >= Nrp*M, "not enough doppler spreading samples");

    hf_gain = 1.0/sqrt(var(spread1)+var(spread2));
    % printf("nsymb: %d lspread1: %d\n", nsymb, length(spread1));
  end
 
  % init timing est states

  tstates.Nsamperframe = Nsamperframe; tstates.M = M; tstates.Ncp = Ncp;
  tstates.verbose = 0; tstates.Fs = Fs;
  delta_t = []; 
  window_width = 11;     % search window_width/2 from current timing instant

  % simulate for each Eb/No point ------------------------------------

  for nn=1:length(EbNodB)
    rand('seed',1);
    randn('seed',1);

    EbNo = bps * (10 .^ (EbNodB(nn)/10));
    variance = 1/(M*EbNo/2);

    % generate tx bits

    tx_bits = rand(1,Nbits) > 0.5;

    % map to symbols in linear array

    if bps == 1
      tx_sym_lin = 2*tx_bits - 1;
    end
    if bps == 2
      for s=1:Nbits/bps
        tx_sym_lin(s) = qpsk_mod(tx_bits(2*(s-1)+1:2*s));
      end
    end

    % place symbols in multi-carrier frame with pilots and boundary carriers

    tx_sym = []; s = 1;
    for f=1:Nframes
      aframe = zeros(Nrowsperframe,Nc+2);
      aframe(1,:) = pilots;
      for r=1:Nrowsperframe
        arowofsymbols = tx_sym_lin(s:s+Nc-1);
        s += Nc;
        aframe(r+1,2:Nc+1) = arowofsymbols;
      end
      tx_sym = [tx_sym; aframe];
    end      
    tx_sym = [tx_sym; pilots];  % final row of pilots
    [nr nc] = size(tx_sym);
    assert(nr == Nrp);

    % OFDM up conversion and upsampling to rate Fs ---------------------

    w = (0:Nc+1)*2*pi*Rs/Fs;
    W = zeros(Nc+2,M);
    for c=1:Nc+2
      W(c,:) = exp(j*w(c)*(0:M-1));
    end

    Nsam = Nrp*(M+Ncp);
    tx = [];

    % OFDM upconvert symbol by symbol so we can add CP

    for r=1:Nrp
      asymbol = tx_sym(r,:) * W/M;
      asymbol_cp = [asymbol(M-Ncp+1:M) asymbol];
      tx = [tx asymbol_cp];
    end
    assert(length(tx) == Nsam);

    % OFDM symbol of pilots is used for coarse timing and freq during acquisition, and fine timing

    rate_fs_pilot_samples = tx(1:M+Ncp);

    % channel simulation ---------------------------------------------------------------

    rx = tx;

    if hf_en
      %rx  =  [zeros(1,path_delay_samples) tx(1:Nsam-path_delay_samples)];
      rx  = hf_gain * tx(1:Nsam) .* spread1(1:Nsam);
      rx  += hf_gain * [zeros(1,path_delay_samples) tx(1:Nsam-path_delay_samples)] .* spread2(1:Nsam);
    end

    rx = rx .* exp(j*woffset*(1:Nsam));

    noise = sqrt(variance)*(0.5*randn(1,Nsam) + j*0.5*randn(1,Nsam));
    rx += noise;

    % some spare samples at end to allow for timing est window

    rx = [rx zeros(1,M)];
    
    % pilot based phase est, we use known tx symbols as pilots ----------
 
    rx_sym = zeros(Nrp, Nc+2);
    phase_est_pilot_log = 10*ones(Nrp,Nc+2);
    phase_est_stripped_log = 10*ones(Nrp,Nc+2);
    phase_est_log = 10*ones(Nrp,Nc+2);
    timing_est_log = [];
    timing_est = Ncp/2;

    for r=1:Ns:Nrp-Ns

      if timing_en

        % update timing every frame

        if (mod(r-1,Ns) == 0) && (r != 1) && (r != Nrp)
          %st = (r-1)*(M+Ncp) - ceil(window_width/2);
          st = (r-1)*(M+Ncp) + 1 - floor(window_width/2) + (timing_est-1);
          en = st + Nsamperframe-1 + length(rate_fs_pilot_samples) + window_width-1;
          ft_est = coarse_sync(tstates, rx(st:en), rate_fs_pilot_samples);
          timing_est += ft_est - ceil(window_width/2);
          %if verbose
          %  printf("ft_est: %d timing_est %d\n", ft_est, timing_est);
          %end
          delta_t = [delta_t ft_est - ceil(window_width/2)];
        end
      end

      timing_est_log = [timing_est_log timing_est];

      % down convert at current timing instant

      % previous pilot

      if r > Ns+1
        rr = r-Ns;
        st = (rr-1)*(M+Ncp) + 1 + timing_est; en = st + M - 1;
        for c=1:Nc+2
          acarrier = rx(st:en) .* conj(W(c,:));
          rx_sym(rr,c) = sum(acarrier);
        end
      end

      % pilot - this frame - pilot

      for rr=r:r+Ns
 
        st = (rr-1)*(M+Ncp) + 1 + timing_est; en = st + M - 1;
        for c=1:Nc+2
          acarrier = rx(st:en) .* conj(W(c,:));
          rx_sym(rr,c) = sum(acarrier);
        end
      end

      % next pilot

      if r < Nrp - 2*Ns
        rr = r+2*Ns;
        st = (rr-1)*(M+Ncp) + 1 + timing_est; en = st + M - 1;
        for c=1:Nc+2
          acarrier = rx(st:en) .* conj(W(c,:));
          rx_sym(rr,c) = sum(acarrier);
        end
      end

      % OK - now estimate and correct phase 

      for c=2:Nc+1

        % estimate phase using average of 6 pilots in a rect 2D window centred
        % on this carrier
        % PPP
        % DDD
        % DDD
        % PPP
          
        cr = c-1:c+1;
        aphase_est_pilot_rect = sum(rx_sym(r,cr)*tx_sym(r,cr)') + sum(rx_sym(r+Ns,cr)*tx_sym(r+Ns,cr)');

        % use next step of pilots in past and future

        if r > Ns+1
          aphase_est_pilot_rect += sum(rx_sym(r-Ns,cr)*tx_sym(r-Ns,cr)');
        end
        if r < Nrp - 2*Ns
          aphase_est_pilot_rect += sum(rx_sym(r+2*Ns,cr)*tx_sym(r+2*Ns,cr)');
        end

        % correct phase offset using phase estimate

        for rr=r+1:r+Ns-1
          aphase_est_pilot = angle(aphase_est_pilot_rect);
          phase_est_pilot_log(rr,c) = aphase_est_pilot;
          rx_corr(rr,c) = rx_sym(rr,c) * exp(-j*aphase_est_pilot);
        end

      end % c=2:Nc+1
    end % r=1:Ns:Nrp-Ns


    % remove pilots to give us just data symbols and demodulate

    rx_bits = []; rx_np = [];
    for r=1:Nrp
      if mod(r-1,Ns) != 0
        arowofsymbols = rx_corr(r,2:Nc+1);
        rx_np = [rx_np arowofsymbols];
        if bps == 1
          arowofbits = real(arowofsymbols) > 0;
        end
        if bps == 2
          arowofbits = zeros(1,Nc);
          for c=1:Nc
            arowofbits((c-1)*2+1:c*2) = qpsk_demod(arowofsymbols(c));
          end
        end
        rx_bits = [rx_bits arowofbits];
      end
    end
    assert(length(rx_bits) == Nbits);


    % calculate BER stats as a block, after pilots extracted

    errors = xor(tx_bits, rx_bits);
    Nerrs = sum(errors);

    printf("EbNodB: %3.2f BER: %4.3f Nbits: %d Nerrs: %d\n", EbNodB(nn), Nerrs/Nbits, Nbits, Nerrs);

    if verbose
      figure(1)
      plot(real(tx))
      figure(2)
      Tx = abs(fft(tx.*hanning(Nsam)'));
      Tx_dB = 20*log10(Tx);
      dF = Fs/Nsam;
      plot((1:Nsam)*dF, Tx_dB);
      mx = max(Tx_dB);
      axis([0 Fs/2 mx-60 mx])
     
      figure(3); clf; 
      plot(rx_np,'+');
      axis([-2 2 -2 2]);

      
      if hf_en
        figure(4); clf; 
        subplot(211)
        plot(abs(spread1(1:Nsam)));
        %hold on; plot(abs(spread2(1:Nsam)),'g'); hold off;
        subplot(212)
        plot(angle(spread1(1:Nsam)));
      end
      

      figure(5); clf;
      plot(phase_est_log(:,2:Nc+1),'+', 'markersize', 10); 
      hold on; 
      plot(phase_est_pilot_log(:,2:Nc+1),'g+', 'markersize', 5); 

#{
      % todo, work out a way to plot rate Fs hf model phase
      if sim_in.hf_en
        plot(angle(hf_model(:,2:Nc+1)));
      end
#}

      axis([1 Nrp -pi pi]);  

      figure(6); clf;
      subplot(211)
      stem(delta_t)
      subplot(212)
      plot(timing_est_log);
    end

    sim_out.ber(nn) = sum(Nerrs)/Nbits; 
    sim_out.pilot_overhead = 10*log10(Ns/(Ns-1));
    sim_out.M = M; sim_out.Fs = Fs; sim_out.Ncp = Ncp;
    sim_out.Nrowsperframe = Nrowsperframe; sim_out.Nsamperframe = Nsamperframe;
  end
endfunction


function run_single
  Ts = 0.016; 
  sim_in.Tcp = 0.002; 
  sim_in.Rs = 1/Ts; sim_in.bps = 2; sim_in.Nc = 16; sim_in.Ns = 8;

  %sim_in.Nsec = 3*(sim_in.Ns+1)/sim_in.Rs;  % one frame
  sim_in.Nsec = 30;

  sim_in.EbNodB = 6;
  sim_in.verbose = 1;
  sim_in.hf_en = 1;
  sim_in.foff_hz = 0;
  sim_in.timing_en = 0;

  run_sim(sim_in);
end


% Plot BER against Eb/No curves for AWGN and HF

% Target operating point Eb/No for HF is 6dB, as this is where our rate 1/2
% LDPC code gives good results (10% PER, 1% BER).  However this means
% the Eb/No at the input is 10*log(1/2) or 3dB less, so we need to
% make sure phase est works at Eb/No = 6 - 3 = 3dB
%
% For AWGN target is 2dB so -1dB op point.

function run_curves
  Ts = 0.010;
  sim_in.Rs = 1/Ts;
  sim_in.Tcp = 0.002; 
  sim_in.bps = 2; sim_in.Ns = 8; sim_in.Nc = 8; sim_in.verbose = 0;
  sim_in.foff_hz = 0;

  pilot_overhead = (sim_in.Ns-1)/sim_in.Ns;
  cp_overhead = Ts/(Ts+sim_in.Tcp);
  overhead_dB = -10*log10(pilot_overhead*cp_overhead);

  sim_in.hf_en = 0;
  sim_in.Nsec = 30;
  sim_in.EbNodB = -3:5;
  awgn_EbNodB = sim_in.EbNodB;

  awgn_theory = 0.5*erfc(sqrt(10.^(sim_in.EbNodB/10)));
  awgn = run_sim(sim_in);

  sim_in.hf_en = 1;
  sim_in.Nsec = 120;
  sim_in.EbNodB = 1:8;

  EbNoLin = 10.^(sim_in.EbNodB/10);
  hf_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));

  hf = run_sim(sim_in);

  figure(4); clf;
  semilogy(awgn_EbNodB, awgn_theory,'b+-;AWGN theory;');
  hold on;
  semilogy(sim_in.EbNodB, hf_theory,'b+-;HF theory;');
  semilogy(awgn_EbNodB+overhead_dB, awgn_theory,'g+-;AWGN lower bound with pilot + CP overhead;');
  semilogy(sim_in.EbNodB+overhead_dB, hf_theory,'g+-;HF lower bound with pilot + CP overhead;');
  semilogy(awgn_EbNodB+overhead_dB, awgn.ber,'r+-;AWGN sim;');
  semilogy(sim_in.EbNodB+overhead_dB, hf.ber,'r+-;HF sim;');
  hold off;
  axis([-3 8 1E-2 2E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
end


% Run an acquisition test, returning vectors of estimation errors

function [delta_t delta_foff] = acquisition_test(Ntests=10, EbNodB=100, foff_hz=0, hf_en=0, fine_en)

  % generate test signal at a given Eb/No and frequency offset

  Ts = 0.016; 
  sim_in.Tcp = 0.002; 
  sim_in.Rs = 1/Ts; sim_in.bps = 2; sim_in.Nc = 16; sim_in.Ns = 8;

  sim_in.Nsec = Ntests*(sim_in.Ns+1)/sim_in.Rs;

  sim_in.EbNodB = EbNodB;
  sim_in.verbose = 0;
  sim_in.hf_en = hf_en;
  sim_in.foff_hz = foff_hz;

  [sim_out rate_fs_pilot_samples rx] = run_sim(sim_in);

  % set up acquistion 

  Nsamperframe = states.Nsamperframe = sim_out.Nsamperframe;
  states.M = sim_out.M; states.Ncp = sim_out.Ncp;
  states.verbose = 0;
  states.Fs = sim_out.Fs;

  % test fine or acquisition over test signal
  #{
    fine: - start with coarse timing instant
          - on each frame est timing a few samples about that point
          - update timing instant

    corr: - where is best plcase to sample
          - just before end of symbol?
          - how long should sequence be?
          - add extra?
          - aim for last possible moment?
          - man I hope IL isn't too big.....
  #}

  delta_t = []; delta_t = [];  delta_foff = [];

  if fine_en

    window_width = 5;                   % search +/-2 samples from current timing instant
    timing_instant = Nsamperframe+1;    % start at correct instant for AWGN
                                        % start at second frame so we can search -2 ... +2

    while timing_instant < (length(rx) - (Nsamperframe + length(rate_fs_pilot_samples) + window_width))
      st = timing_instant - ceil(window_width/2); 
      en = st + Nsamperframe-1 + length(rate_fs_pilot_samples) + window_width-1;
      [ft_est foff_est] = coarse_sync(states, rx(st:en), rate_fs_pilot_samples);
      printf("ft_est: %d timing_instant %d %d\n", ft_est, timing_instant, mod(timing_instant, Nsamperframe));
      timing_instant += Nsamperframe + ft_est - ceil(window_width/2);
      delta_t = [delta_ft ft_est - ceil(window_width/2)];
    end
  else
    % for coarse simulation we just use contant window shifts

    st = 0.5*Nsamperframe; 
    en = 2.5*Nsamperframe - 1;
    ct_target = Nsamperframe/2;

    for w=1:Nsamperframe:length(rx)-3*Nsamperframe
      %st = w+0.5*Nsamperframe; en = st+2*Nsamperframe-1;
      %[ct_est foff_est] = coarse_sync(states, rx(st:en), rate_fs_pilot_samples);
      [ct_est foff_est] = coarse_sync(states, rx(w+st:w+en), rate_fs_pilot_samples);
      printf("ct_est: %4d foff_est: %3.1f\n", ct_est, foff_est);

      % valid coarse timing ests are modulo Nsamperframe

      delta_t = [delta_ct ct_est-ct_target];
      delta_foff = [delta_foff (foff_est-foff_hz)];
    end
  end

endfunction


% Run some tests for various acquisition conditions. Probability of
% acquistion is what matters, e.g. if it's 50% we can expect sync
% within 2 frames
%                 P(t)/P(f)  P(t)/P(f)
%          Eb/No  AWGN       HF
% +/- 25Hz -1/3   1.0/0.3    0.96/0.3 
% +/-  5Hz -1/3   1.0/0.347  0.96/0.55 
% +/- 25Hz  10/10 1.00/0.92  0.99/0.77

function acquisition_histograms(fine_en = 0)
  Fs = 8000;
  Ntests = 100;

  % allowable tolerance for acquistion

  ftol_hz = 2.0;
  ttol_samples = 0.002*Fs;

  % AWGN channel operating point

  [dct dfoff] = acquisition_test(Ntests, -1, 25, 0, fine_en);

  % Probability of acquistion is what matters, e.g. if it's 50% we can
  % expect sync within 2 frames

  printf("AWGN P(time offset acq) = %3.2f\n", length(find (abs(dct) < ttol_samples))/length(dct));
  if fine_en == 0
    printf("AWGN P(freq offset acq) = %3.2f\n", length(find (abs(dfoff) < ftol_hz))/length(dfoff));
  end

  figure(1)
  hist(dct(find (abs(dct) < ttol_samples)))
  if fine_en == 0
    figure(2)
    hist(dfoff)
  end

  % HF channel operating point

  [dct dfoff] = acquisition_test(Ntests, 3, 25, 1, fine_en);

  printf("HF P(time offset acq) = %3.2f\n", length(find (abs(dct) < ttol_samples))/length(dct));
  if fine_en == 0
    printf("HF P(freq offset acq) = %3.2f\n", length(find (abs(dfoff) < ftol_hz))/length(dfoff));
  end

  figure(3)
  hist(dct(find (abs(dct) < ttol_samples)))
  if fine_en == 0
    figure(4)
    hist(dfoff)
  end

endfunction



% choose simulation to run here -------------------------------------------------------

format;
more off;

run_single
%run_curves
%acquisition_histograms(1)

