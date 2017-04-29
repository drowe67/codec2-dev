% ofdm_fs.m
% David Rowe Mar 2017
%
% Rate Fs BPSK/QPSK OFDM simulation, rate Fs verison of ofdm_rs.m with
% OFDM based up and down conversion.

#{
 TODO:
   [X] strip back experimental stuff to just features we need
   [X] ZOH/integrator
   [X] OFDM up and down conversion
   [X] rate Fs HF model and HF results
   [X] add QPSK
   [X] add CP
   [X] fine timing estimator and sample clock offset tracking
   [X] acquisition coarse timing & freq offset estimation
   [ ] adjust waveform parameters for real world
   [ ] Nsec run time take into account CP
   [ ] plot phase ests
   [ ] handle border carriers
       [ ] start with phantom carriers 
           + but unhappy with 1800Hz bandwidth
       [ ] also try interpolation or just single row
   [ ] compute SNR and PAPR
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


% Correlates the OFDM pilot symbol samples with a window of received
% samples to determine the most likely timing offset.  Combines two
% frames pilots so we need at least Nsamperframe+M+Ncp samples in rx.
% Also determines frequency offset at maximimum correlation.  Can be
% used for acquisition (coarse timing a freq offset), and fine timing

function [t_est foff_est] = coarse_sync(states, rx, rate_fs_pilot_samples)
    Nsamperframe = states.Nsamperframe; Fs = states.Fs;
    Npsam = length(rate_fs_pilot_samples);
    verbose  = states.verbose;

    Ncorr = length(rx) - (Nsamperframe+Npsam) + 1;
    assert(Ncorr > 0);
    corr = zeros(1,Ncorr);
    for i=1:Ncorr
      corr(i)   = abs(rx(i:i+Npsam-1) * rate_fs_pilot_samples');
      corr(i)  += abs(rx(i+Nsamperframe:i+Nsamperframe+Npsam-1) * rate_fs_pilot_samples');
    end

    [mx t_est] = max(abs(corr));

    C  = abs(fft(rx(t_est:t_est+Npsam-1) .* conj(rate_fs_pilot_samples), Fs));
    C += abs(fft(rx(t_est+Nsamperframe:t_est+Nsamperframe+Npsam-1) .* conj(rate_fs_pilot_samples), Fs));

    fmax = 30;
    [mx_pos foff_est_pos] = max(C(1:fmax));
    [mx_neg foff_est_neg] = max(C(Fs-fmax+1:Fs));

    if mx_pos > mx_neg
      foff_est = foff_est_pos - 1;
    else
      foff_est = foff_est_neg - fmax - 1;
    end

    if verbose > 1
      %printf("t_est: %d\n", t_est);
      figure(7); clf;
      plot(abs(corr))
      figure(8)
      plot(C)
      axis([0 200 0 0.4])
    end
endfunction


#{
  Frame has Ns-1 data symbols between pilots, e.g. for Ns=3: 
  
    PPP
    DDD
    DDD
    PPP
#}

function states = ofdm_init(bps, Rs, Tcp, Ns, Nc)
  states.Fs = 8000;
  states.bps = bps;
  states.Rs = Rs;
  states.Tcp = Tcp;
  states.Ns = Ns;       % step size for pilots
  states.Nc = Nc;       % Number of cols, aka number of carriers
  states.M  = states.Fs/Rs;
  states.Ncp = Tcp*states.Fs;
  states.Nbitsperframe = (Ns-1)*Nc*bps;
  states.Nrowsperframe = states.Nbitsperframe/(Nc*bps);
  states.Nsamperframe =  (states.Nrowsperframe+1)*(states.M+states.Ncp);

  % generate same pilots each time

  rand('seed',1);
  states.pilots = 1 - 2*(rand(1,Nc+2) > 0.5);

  % carrier tables for up and down conversion

  w = (0:Nc+1)*2*pi*Rs/states.Fs;
  W = zeros(Nc+2,states.M);
  for c=1:Nc+2
    W(c,:) = exp(j*w(c)*(0:states.M-1));
  end
  states.w = w;
  states.W = W;

  % fine timing search +/- window_width/2 from current timing instant

  states.window_width = 11; 
 
  % Receive buffer: D P DDD P DDD P DDD P D
  %                         ^
  % also see ofdm_demod() ...

  states.Nrxbuf = 3*states.Nsamperframe+states.M+states.Ncp + 2*(states.M + states.Ncp);
  states.rxbuf = zeros(1, states.Nrxbuf);
 
  % default settings on a bunch of options and states

  states.verbose = 0;
  states.timing_en = 1;
  states.foff_est_en = 1;
  states.phase_est_en = 1;

  states.foff_est_gain = 0.1;
  states.foff_est_hz = 0;
  states.sample_point = states.timing_est = 1;
  states.nin = states.Nsamperframe;

endfunction


% ---------------------------
% Modulates one frame of bits
% ---------------------------

function tx = ofdm_mod(states, tx_bits)
  ofdm_load_const;

  assert(length(tx_bits) == Nbitsperframe);

  % map to symbols in linear array

  if bps == 1
    tx_sym_lin = 2*tx_bits - 1;
  end
  if bps == 2
    for s=1:Nbitsperframe/bps
      tx_sym_lin(s) = qpsk_mod(tx_bits(2*(s-1)+1:2*s));
    end
  end

  % place symbols in multi-carrier frame with pilots and boundary carriers

  tx_sym = []; s = 1;
  aframe = zeros(Ns,Nc+2);
  aframe(1,:) = pilots;
  for r=1:Nrowsperframe
    arowofsymbols = tx_sym_lin(s:s+Nc-1);
    s += Nc;
    aframe(r+1,2:Nc+1) = arowofsymbols;
  end
  tx_sym = [tx_sym; aframe];

  % OFDM upconvert symbol by symbol so we can add CP

  tx = [];
  for r=1:Ns
    asymbol = tx_sym(r,:) * W/M;
    asymbol_cp = [asymbol(M-Ncp+1:M) asymbol];
    tx = [tx asymbol_cp];
  end
endfunction


% -----------------------------
% Demodulates one frame of bits
% -----------------------------

#{ 

  For phase estimation we need to maintain buffer of 3 frames plus
  one pilot, so we have 4 pilots total. '^' is the start of current
  frame that we are demodulating.
           
  P DDD P DDD P DDD P
        ^
    
  Then add one symbol either side to account for movement in
  sampling instant due to sample clock differences:

  D P DDD P DDD P DDD P D
          ^
#}

function [rx_bits states aphase_est_pilot_log rx_np] = ofdm_demod(states, rxbuf_in)
  ofdm_load_const;

  % extra states that are st up at run time rather than init time

  timing_est = states.timing_est;
  timing_en = states.timing_en;
  foff_est_hz = states.foff_est_hz;
  foff_est_gain = states.foff_est_gain;
  foff_est_en = states.foff_est_en;
  sample_point = states.sample_point;
  rate_fs_pilot_samples = states.rate_fs_pilot_samples;
  verbose = states.verbose;
  phase_est_en = states.phase_est_en;

  % insert latest input samples into rxbuf

  rxbuf(1:Nrxbuf-states.nin) = rxbuf(states.nin+1:Nrxbuf);
  rxbuf(Nrxbuf-states.nin+1:Nrxbuf) = rxbuf_in;

  woff_est = 2*pi*foff_est_hz/Fs;

  % update timing estimate --------------------------------------------------

  delta_t = sample_point = 0;
  if timing_en
    % update timing at start of every frame

    st = M+Ncp + Nsamperframe + 1 - floor(window_width/2) + (timing_est-1);
    en = st + Nsamperframe-1 + M+Ncp + window_width-1;
          
    ft_est = coarse_sync(states, rxbuf(st:en) .* exp(-j*woff_est*(st:en)), rate_fs_pilot_samples);
    timing_est += ft_est - ceil(window_width/2);

    if verbose > 1
      printf("  ft_est: %2d timing_est: %2d sample_point: %2d\n", ft_est, timing_est, sample_point);
    end

    % Black magic to keep sample_point inside cyclic prefix.  Or something like that.

    delta_t = ft_est - ceil(window_width/2);
    sample_point = max(timing_est+Ncp/4, sample_point);
    sample_point = min(timing_est+Ncp, sample_point);
  end

  % down convert at current timing instant----------------------------------

    % todo: this cld be more efficent, as pilot r becomes r-Ns on next frame

  rx_sym = zeros(1+Ns+1+1, Nc+2);

  % previous pilot

  st = M+Ncp + Nsamperframe + (-Ns)*(M+Ncp) + 1 + sample_point; en = st + M - 1;

  for c=1:Nc+2
    acarrier = rxbuf(st:en) .* exp(-j*woff_est*(st:en)) .* conj(W(c,:));
    rx_sym(1,c) = sum(acarrier);
  end

  % pilot - this frame - pilot

  for rr=1:Ns+1 
    st = M+Ncp + Nsamperframe + (rr-1)*(M+Ncp) + 1 + sample_point; en = st + M - 1;
    for c=1:Nc+2
      acarrier = rxbuf(st:en) .* exp(-j*woff_est*(st:en)) .* conj(W(c,:));
      rx_sym(rr+1,c) = sum(acarrier);
    end
  end

  % next pilot

  st = M+Ncp + Nsamperframe + (2*Ns)*(M+Ncp) + 1 + sample_point; en = st + M - 1;
  for c=1:Nc+2
    acarrier = rxbuf(st:en) .* exp(-j*woff_est*(st:en)) .* conj(W(c,:));
    rx_sym(Ns+3,c) = sum(acarrier);
  end
      
  % est freq err based on all carriers ------------------------------------
      
  if foff_est_en
    freq_err_rect = sum(rx_sym(2,:))' * sum(rx_sym(2+Ns,:));
    freq_err_hz = angle(freq_err_rect)*Rs/(2*pi*Ns);
    foff_est_hz += foff_est_gain*freq_err_hz;
  end

  % OK - now estimate and correct phase  ----------------------------------

  aphase_est_pilot = 10*ones(1,Nc+2);
  for c=2:Nc+1

    % estimate phase using average of 6 pilots in a rect 2D window centred
    % on this carrier
    % PPP
    % DDD
    % DDD
    % PPP
          
    cr = c-1:c+1;
    aphase_est_pilot_rect = sum(rx_sym(2,cr)*pilots(cr)') + sum(rx_sym(2+Ns,cr)*pilots(cr)');

    % use next step of pilots in past and future

    aphase_est_pilot_rect += sum(rx_sym(1,cr)*pilots(cr)');
    aphase_est_pilot_rect += sum(rx_sym(2+Ns+1,cr)*pilots(cr)');

    aphase_est_pilot(c) = angle(aphase_est_pilot_rect);
  end

  % correct phase offset using phase estimate, and demodulate
  % bits, separate loop as it runs across cols (carriers) to get
  % frame bit ordering correct

  aphase_est_pilot_log = [];
  rx_bits = []; rx_np = [];
  for rr=1:Ns-1
    for c=2:Nc+1
      if phase_est_en
        rx_corr = rx_sym(rr+2,c) * exp(-j*aphase_est_pilot(c));
      else
        rx_corr = rx_sym(rr+2,c);
      end
      rx_np = [rx_np rx_corr];
      if bps == 1
        abit = real(rx_corr) > 0;
      end
      if bps == 2
        abit = qpsk_demod(rx_corr);
      end
      rx_bits = [rx_bits abit];
    end % c=2:Nc+1
    aphase_est_pilot_log = [aphase_est_pilot_log; aphase_est_pilot];
  end 

  % Adjust nin to take care of sample clock offset

  nin = Nsamperframe;
  if timing_en
    thresh = (M+Ncp)/8;
    tshift = (M+Ncp)/4;
    if timing_est > thresh
      nin = Nsamperframe+tshift;
      timing_est -= tshift;
      sample_point -= tshift;
    end
    if timing_est < -thresh
      nin = Nsamperframe-tshift;
      timing_est += tshift;
      sample_point += tshift;
    end
  end

  states.rxbuf = rxbuf;
  states.nin = nin;
  states.timing_est = timing_est;
  states.sample_point = sample_point;
  states.delta_t = delta_t;
  states.foff_est_hz = foff_est_hz;
endfunction


#{
  TODO: 
    [ ] Some states need warm rest at the start of each simulation point
    [ ] rate_fs_pilot_samples generated in init
    [ ] move buffer shift into demod
    [ ] way to simulate aquisition and demod
    [ ] testframe based 
#}

function [sim_out rate_fs_pilot_samples rx] = run_sim(sim_in)

  % set up core modem constants

  states = ofdm_init(sim_in.bps, sim_in.Rs, sim_in.Tcp, sim_in.Ns, sim_in.Nc);
  ofdm_load_const;
  
  % simulation parameters and flags

  woffset = 2*pi*sim_in.foff_hz/Fs;
  EbNodB  = sim_in.EbNodB;
  verbose = states.verbose = sim_in.verbose;
  hf_en   = sim_in.hf_en;
  timing_en = states.timing_en = sim_in.timing_en;
  states.foff_est_en = foff_est_en = sim_in.foff_est_en;
  states.phase_est_en = phase_est_en = sim_in.phase_est_en;

  if verbose
    printf("Rs:..........: %4.2f\n", Rs);
    printf("M:...........: %d\n", M);
    printf("Ncp:.........: %d\n", Ncp);
    printf("bps:.........: %d\n", bps);
    printf("Nbitsperframe: %d\n", Nbitsperframe);
    printf("Nrowsperframe: %d\n", Nrowsperframe);
    printf("Nsamperframe.: %d\n", Nsamperframe);
  end

  % Important to define run time in seconds so HF model will evolve the same way
  % for different pilot insertion rates.  So lets work backwards from approx
  % seconds in run to get Nbits, the total number of payload data bits

  Nrows = sim_in.Nsec*Rs;
  Nframes = floor((Nrows-1)/Ns);
  Nbits = Nframes * Nbitsperframe;    % number of payload data bits

  Nr = Nbits/(Nc*bps);                % Number of data rows to get Nbits total

  % double check if Nbits fit neatly into carriers

  assert(Nbits/(Nc*bps) == floor(Nbits/(Nc*bps)), "Nbits/(Nc*bps) must be an integer");

  Nrp = Nr + Nframes + 1;  % number of rows once pilots inserted
                           % extra row of pilots at end

  if verbose
    printf("Nc.....: %d\n", Nc);
    printf("Ns.....: %d (step size for pilots, Ns-1 data symbols between pilots)\n", Ns);
    printf("Nr.....: %d\n", Nr);
    printf("Nbits..: %d\n", Nbits);
    printf("Nframes: %d\n", Nframes);
    printf("Nrp....: %d (number of rows including pilots)\n", Nrp);
  end

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
 
  % ------------------------------------------------------------------
  % simulate for each Eb/No point 
  % ------------------------------------------------------------------

  for nn=1:length(EbNodB)
    rand('seed',1);
    randn('seed',1);

    EbNo = bps * (10 .^ (EbNodB(nn)/10));
    variance = 1/(M*EbNo/2);

    Nsam = Nrp*(M+Ncp);

    % generate tx bits and run OFDM modulator

    tx_bits = rand(1,Nbits) > 0.5;

    tx = [];
    for f=1:Nframes
      tx = [tx ofdm_mod(states, tx_bits((f-1)*Nbitsperframe+1:f*Nbitsperframe))];
    end

    % add extra row of pilots at end, to allow one frame simulations,
    % useful for development

    st = Nsamperframe*(Nframes-1)+1; en = st+Ncp+M-1;
    tx = [tx tx(st:en)];
    assert(length(tx) == Nsam);

    % OFDM symbol of pilots is used for coarse timing and freq during acquisition, and fine timing
    % TODO: put this in init code

    states.rate_fs_pilot_samples = tx(1:M+Ncp);

    % channel simulation ---------------------------------------------------------------
    
    if isfield(sim_in, "sample_clock_offset_ppm") 
      % todo: this only works for large ppm like 500, runs out of memory
      %       for small ppm

      if sim_in.sample_clock_offset_ppm
        timebase = floor(abs(1E6/sim_in.sample_clock_offset_ppm));
        if sim_in.sample_clock_offset_ppm > 0
          tx = resample(tx, timebase+1, timebase);
        else
          tx = resample(tx, timebase, timebase+1);
        end

        % make sure length is correct for rest of simulation

        tx = [tx zeros(1,Nsam-length(tx))];
        tx = tx(1:Nsam);
      end
    end

    rx = tx;

    if hf_en
      %rx  =  [zeros(1,path_delay_samples) tx(1:Nsam-path_delay_samples)];
      rx  = hf_gain * tx(1:Nsam) .* spread1(1:Nsam);
      rx  += hf_gain * [zeros(1,path_delay_samples) tx(1:Nsam-path_delay_samples)] .* spread2(1:Nsam);
    end

    rx = rx .* exp(j*woffset*(1:Nsam));

    noise = sqrt(variance)*(0.5*randn(1,Nsam) + j*0.5*randn(1,Nsam));
    rx += noise;

    % some spare samples at end to avoid overflow as est windows may poke into the future a bit

    rx = [rx zeros(1,Nsamperframe)];
    
    % bunch of logs

    phase_est_pilot_log = [];
    delta_t_log = []; 
    timing_est_log = [];
    foff_est_hz_log = [];
    Nerrs_log = [];
    rx_bits = []; rx_np = [];

    % reset some states for each EbNo simulation point

    states.sample_point = states.timing_est = 1;
    if timing_en == 0
      states.sample_point = Ncp;
    end
    states.nin = Nsamperframe;
    states.foff_est_hz = 0;    

    % for this simulation we "prime" buffer to allow one frame runs during development

    prx = 1;
    states.rxbuf(M+Ncp+2*Nsamperframe+1:Nrxbuf) = rx(prx:Nsamperframe+2*(M+Ncp));
    prx += Nsamperframe+2*(M+Ncp);

    for f=1:Nframes

      % insert samples at end of buffer, set to zero if no samples
      % available to disable phase estimation on future pilots on last
      % frame of simulation

      lnew = min(Nsam-prx,states.nin);
      rxbuf_in = zeros(1,states.nin);

      if lnew
        rxbuf_in(1:lnew) = rx(prx:prx+lnew-1);
      end
      prx += states.nin;

      [arx_bits states aphase_est_pilot_log arx_np] = ofdm_demod(states, rxbuf_in);

      rx_bits = [rx_bits arx_bits]; rx_np = [rx_np arx_np];
      timing_est_log = [timing_est_log states.timing_est];
      delta_t_log = [delta_t_log states.delta_t];
      foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
      phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];
    end

    assert(length(rx_bits) == Nbits);

    % calculate BER stats as a block, after pilots extracted

    errors = xor(tx_bits, rx_bits);
    Nerrs = sum(errors);
    for f=1:Nframes
      st = (f-1)*Nbitsperframe+1; en = st + Nbitsperframe-1;
      Nerrs_log(f) = sum(xor(tx_bits(st:en), rx_bits(st:en)));
    end

    printf("EbNodB: %3.2f BER: %5.4f Nbits: %d Nerrs: %d\n", EbNodB(nn), Nerrs/Nbits, Nbits, Nerrs);

    if verbose

      figure(1); clf; 
      plot(rx_np,'+');
      axis([-2 2 -2 2]);
      title('Scatter');
      
      figure(2); clf;
      plot(phase_est_pilot_log(:,2:Nc+1),'g+', 'markersize', 5); 
      title('Phase est');
      axis([1 Nrp -pi pi]);  

      figure(3); clf;
      subplot(211)
      stem(delta_t_log)
      title('delta t');
      subplot(212)
      plot(timing_est_log);
      title('timing est');

      figure(4); clf;
      plot(foff_est_hz_log)
      axis([1 max(Nframes,2) -3 3]);
      title('Fine Freq');

      figure(5); clf;
      plot(Nerrs log);

#{
      figure(2)
      Tx = abs(fft(tx(1:Nsam).*hanning(Nsam)'));
      Tx_dB = 20*log10(Tx);
      dF = Fs/Nsam;
      plot((1:Nsam)*dF, Tx_dB);
      mx = max(Tx_dB);
      axis([0 Fs/2 mx-60 mx])
#}
     
#{
      if hf_en
        figure(4); clf; 
        subplot(211)
        plot(abs(spread1(1:Nsam)));
        %hold on; plot(abs(spread2(1:Nsam)),'g'); hold off;
        subplot(212)
        plot(angle(spread1(1:Nsam)));
        title('spread1 amp and phase');
      end
#}
      
#{
      % todo, work out a way to plot rate Fs hf model phase
      if sim_in.hf_en
        plot(angle(hf_model(:,2:Nc+1)));
      end
#}


    end

    sim_out.ber(nn) = sum(Nerrs)/Nbits; 
    sim_out.pilot_overhead = 10*log10(Ns/(Ns-1));
    sim_out.M = M; sim_out.Fs = Fs; sim_out.Ncp = Ncp;
    sim_out.Nrowsperframe = Nrowsperframe; sim_out.Nsamperframe = Nsamperframe;
  end
endfunction


function run_single
  Ts = 0.018; 
  sim_in.Tcp = 0.002; 
  sim_in.Rs = 1/Ts; sim_in.bps = 2; sim_in.Nc = 16; sim_in.Ns = 8;

  %sim_in.Nsec = 5*(sim_in.Ns+1)/sim_in.Rs;  % one frame
  sim_in.Nsec = 30;

  sim_in.EbNodB = 3;
  sim_in.verbose = 1;
  sim_in.hf_en = 1;
  sim_in.foff_hz = 0;
  sim_in.timing_en = 1;
  sim_in.sample_clock_offset_ppm = 0;
  sim_in.foff_est_en = 1;
  sim_in.phase_est_en = 1;

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

function [delta_t delta_foff] = acquisition_test(Ntests=10, EbNodB=100, foff_hz=0, hf_en=0, fine_en=0)

  % generate test signal at a given Eb/No and frequency offset

  Ts = 0.016; 
  sim_in.Tcp = 0.002; 
  sim_in.Rs = 1/Ts; sim_in.bps = 2; sim_in.Nc = 16; sim_in.Ns = 8;

  sim_in.Nsec = Ntests*(sim_in.Ns+1)/sim_in.Rs;

  sim_in.EbNodB = EbNodB;
  sim_in.verbose = 0;
  sim_in.hf_en = hf_en;
  sim_in.foff_hz = foff_hz; sim_in.timing_en = 0;

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

