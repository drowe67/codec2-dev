% ofdm_lib.m
% David Rowe Mar 2017

#{
  Library of functions that implement a BPSK/QPSK OFDM modem.  Rate Fs
  verison of ofdm_rs.m with OFDM based up and down conversion, and all
  those nasty real-world details like fine freq, timing.  
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

function out = freq_shift(in, foff, Fs)
  foff_rect = exp(j*2*pi*foff/Fs);
  foff_phase_rect = exp(j*0);
  
  for r=1:length(in)
    foff_phase_rect *= foff_rect;
    out(r) = in(r)*foff_phase_rect;
  end
endfunction

% Correlates the OFDM pilot symbol samples with a window of received
% samples to determine the most likely timing offset.  Combines two
% frames pilots so we need at least Nsamperframe+M+Ncp samples in rx.
% Also determines frequency offset at maximimum correlation.  Can be
% used for acquisition (coarse timing and freq offset), and fine
% timing

function [t_est foff_est timing_valid timing_mx] = coarse_sync(states, rx, rate_fs_pilot_samples)
    ofdm_load_const;
    Npsam = length(rate_fs_pilot_samples);

    Ncorr = length(rx) - (Nsamperframe+Npsam) + 1;
    assert(Ncorr > 0);
    corr = zeros(1,Ncorr);

    % normalise correlation so we can compare to a threshold across varying input levels

    av_level = 2*sqrt(states.timing_norm*(rx*rx')/length(rx)) + 1E-12;

    % correlate with pilots at start and end of frame to determine timing offset
    
    for i=1:Ncorr
      rx1     = rx(i:i+Npsam-1); rx2 = rx(i+Nsamperframe:i+Nsamperframe+Npsam-1);
      corr_st = rx1 * rate_fs_pilot_samples'; corr_en = rx2 * rate_fs_pilot_samples';
      corr(i) = (abs(corr_st) + abs(corr_en))/av_level;
    end

    [timing_mx t_est] = max(corr);
    timing_valid = timing_mx > timing_mx_thresh;
    
    if verbose > 2
      printf("   mx: %f timing_est: %d timing_valid: %d\n", timing_mx, t_est, timing_valid);
    end
        
    p1 = rx(t_est:t_est+Npsam/2-1) * rate_fs_pilot_samples(1:Npsam/2)';
    p2 = rx(t_est+Npsam/2:t_est+Npsam-1) * rate_fs_pilot_samples(Npsam/2+1:Npsam)';
    p3 = rx(t_est+Nsamperframe:t_est+Nsamperframe+Npsam/2-1) * rate_fs_pilot_samples(1:Npsam/2)';
    p4 = rx(t_est+Nsamperframe+Npsam/2:t_est+Nsamperframe+Npsam-1) * rate_fs_pilot_samples(Npsam/2+1:Npsam)';
   
    Fs1 = Fs/(Npsam/2);
    foff_est = Fs1*angle(conj(p1)*p2 + conj(p3)*p4)/(2*pi);
        
    if verbose > 1
      figure(7); clf;
      plot(abs(corr))

      figure(8)
      subplot(211)
      plot(real(rate_fs_pilot_samples))
      hold on; plot(real(rx(t_est:t_est+Npsam-1)),'g'); hold off
      subplot(212)
      plot(imag(rate_fs_pilot_samples))
      hold on; plot(imag(rx(t_est:t_est+Npsam-1)),'g'); hold off

      figure(9)
      dp = rx(t_est:t_est+Npsam-1) .* conj(rate_fs_pilot_samples);
      subplot(211); plot(real(dp));
      subplot(212); plot(imag(dp));

      figure(10)
      plot(dp,'+')
      
      figure(11)
      plot([p1 p2; p3 p4], '+')
      
    end

endfunction


%-------------------------------------------------------------
% ofdm_init
%-------------------------------------------------------------

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
  states.Ntxtbits = 4;   % reserve 4 bits/frame for auxillary text information
  states.Nuwbits  = (Ns-1)*bps - states.Ntxtbits;
  states.tx_uw = [1 0 0 1 0 1 0 0 1 0];
  assert(length(states.tx_uw) == states.Nuwbits);

  % use this to scale tx output to 16 bit short.  Adjusted by experiment
  % to have same RMS power as FDMDV waveform
  
  states.amp_scale = 2E5*1.1491/1.06;

  % this is used to scale input sto LDPC decoder to make it amplitude indep
  
  states.mean_amp = 0;

  % generate same pilots each time

  rand('seed',1);
  states.pilots = 1 - 2*(rand(1,Nc+2) > 0.5);
  
  % carrier tables for up and down conversion

  fcentre = 1500;
  Nlower = floor((fcentre - Rs*Nc/2)/Rs);
  w = (Nlower:Nlower+Nc+1)*2*pi*Rs/states.Fs;
  W = zeros(Nc+2,states.M);
  for c=1:Nc+2
    W(c,:) = exp(j*w(c)*(0:states.M-1));
  end
  states.w = w;
  states.W = W;

  % fine timing search +/- window_width/2 from current timing instant

  states.ftwindow_width = 11; 
 
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

  states.foff_est_gain = 0.05;
  states.foff_est_hz = 0;
  states.sample_point = states.timing_est = 1;
  states.nin = states.Nsamperframe;
  states.timing_valid = 0;
  states.timing_mx = 0;
  states.coarse_foff_est_hz = 0;
  
  % generate OFDM pilot symbol, used for timing and freq offset est

  rate_fs_pilot_samples = states.pilots * W/states.M;
  states.rate_fs_pilot_samples = [rate_fs_pilot_samples(states.M-states.Ncp+1:states.M) rate_fs_pilot_samples];

  % pre-compute a constant used in coarse_sync()

  Npsam = length(states.rate_fs_pilot_samples);
  states.timing_norm = Npsam*(states.rate_fs_pilot_samples * states.rate_fs_pilot_samples');
  % printf("timing_norm: %f\n", states.timing_norm)

  % sync state machine

  states.sync_state = states.last_sync_state = 'search';
  states.uw_errors = 0;
  states.sync_counter = 0;
  states.frame_count = 0;
  states.sync_start = 0;
  states.sync_end = 0;
  states.sync_state_interleaver = 'search';
  states.frame_count_interleaver = 0;
   
  % LDPC code is optionally enabled

  states.rate = 1.0;
  states.ldpc_en = 0;

  % init some output states for logging
  
  states.rx_sym = zeros(1+Ns+1+1, Nc+2);

  % Es/No (SNR) est states
  
  states.noise_var = 0;
  states.sig_var = 0;
  
endfunction



% --------------------------------------
% ofdm_mod - modulates one frame of bits
% --------------------------------------

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

  tx = ofdm_txframe(states, tx_sym_lin);
endfunction

% -----------------------------------------
% ofdm_tx - modulates one frame of symbols
% ----------------------------------------

#{ 
   Each carrier amplitude is 1/M.  There are two edge carriers that
   are just tx-ed for pilots plus plus Nc continuous carriers. So
   power is:

     p = 2/(Ns*(M*M)) + Nc/(M*M)

   e.g. Ns=8, Nc=16, M=144 
   
     p = 2/(8*(144*144)) + 16/(144*144) = 7.84-04

#}

function tx = ofdm_txframe(states, tx_sym_lin)
  ofdm_load_const;
  assert(length(tx_sym_lin) == Nbitsperframe/bps);

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


% ----------------------------------------------------------------------------------
% ofdm_sync_search - attempts to find coarse sync parameters for modem initial sync
% ----------------------------------------------------------------------------------

function [timing_valid states] = ofdm_sync_search(states, rxbuf_in)
  ofdm_load_const;

  % insert latest input samples into rxbuf so it is primed for when we have to call ofdm_demod()

  states.rxbuf(1:Nrxbuf-states.nin) = states.rxbuf(states.nin+1:Nrxbuf);
  states.rxbuf(Nrxbuf-states.nin+1:Nrxbuf) = rxbuf_in;

  % Attempt coarse timing estimate (i.e. detect start of frame)

  st = M+Ncp + Nsamperframe + 1; en = st + 2*Nsamperframe; 
  [ct_est foff_est timing_valid timing_mx] = coarse_sync(states, states.rxbuf(st:en), states.rate_fs_pilot_samples);
  if verbose
    printf("  ct_est: %d mx: %3.2f coarse_foff: %4.1f\n", ct_est, timing_mx, foff_est);
  end

  if timing_valid
    % potential candidate found ....

    % calculate number of samples we need on next buffer to get into sync

    states.nin = Nsamperframe + ct_est - 1;

    % reset modem states

    states.sample_point = states.timing_est = 1;
    states.foff_est_hz = foff_est;
  else
    states.nin = Nsamperframe;
  end
  
  states.timing_valid = timing_valid;
  states.timing_mx = timing_mx;
  states.coarse_foff_est_hz = foff_est;
endfunction


% ------------------------------------------
% ofdm_demod - Demodulates one frame of bits
% ------------------------------------------

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

function [rx_bits states aphase_est_pilot_log rx_np rx_amp] = ofdm_demod(states, rxbuf_in)
  ofdm_load_const;

  % insert latest input samples into rxbuf

  rxbuf(1:Nrxbuf-states.nin) = rxbuf(states.nin+1:Nrxbuf);
  rxbuf(Nrxbuf-states.nin+1:Nrxbuf) = rxbuf_in;

  % get latest freq offset estimate

  woff_est = 2*pi*foff_est_hz/Fs;

  % update timing estimate --------------------------------------------------

  delta_t = coarse_foff_est_hz = timing_valid = timing_mx = 0;
  if timing_en
    % update timing at start of every frame

    st = M+Ncp + Nsamperframe + 1 - floor(ftwindow_width/2) + (timing_est-1);
    en = st + Nsamperframe-1 + M+Ncp + ftwindow_width-1;
          
    [ft_est coarse_foff_est_hz timing_valid timing_mx] = coarse_sync(states, rxbuf(st:en) .* exp(-j*woff_est*(st:en)), rate_fs_pilot_samples);

    if timing_valid
      timing_est = timing_est + ft_est - ceil(ftwindow_width/2);

      % Black magic to keep sample_point inside cyclic prefix.  Or something like that.

      delta_t = ft_est - ceil(ftwindow_width/2);
      sample_point = max(timing_est+Ncp/4, sample_point);
      sample_point = min(timing_est+Ncp, sample_point);
    end
    
    if verbose > 1
      printf("  ft_est: %2d mx: %3.2f coarse_foff: %4.1f foff: %4.1f\n", ft_est, timing_mx, coarse_foff_est_hz, foff_est_hz);
    end

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

    % prevent instability in atan(im/re) when real part near 0 

    freq_err_rect += 1E-6;

    %printf("freq_err_rect: %f %f angle: %f\n", real(freq_err_rect), imag(freq_err_rect), angle(freq_err_rect));
    freq_err_hz = angle(freq_err_rect)*Rs/(2*pi*Ns);
    foff_est_hz = foff_est_hz + foff_est_gain*freq_err_hz;
  end

  % OK - now estimate and correct phase  ----------------------------------

  aphase_est_pilot = 10*ones(1,Nc+2);
  aamp_est_pilot = zeros(1,Nc+2);
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

    % amplitude is estimated over 12 pilot symbols, so find average

    aamp_est_pilot(c) = abs(aphase_est_pilot_rect/12);
  end
 
  % correct phase offset using phase estimate, and demodulate
  % bits, separate loop as it runs across cols (carriers) to get
  % frame bit ordering correct

  aphase_est_pilot_log = [];
  rx_bits = []; rx_np = []; rx_amp = [];
  for rr=1:Ns-1
    for c=2:Nc+1
      if phase_est_en
        rx_corr = rx_sym(rr+2,c) * exp(-j*aphase_est_pilot(c));
      else
        rx_corr = rx_sym(rr+2,c);
      end
      rx_np = [rx_np rx_corr];
      rx_amp = [rx_amp aamp_est_pilot(c)];
      if bps == 1
        abit = real(rx_corr) > 0;
      end
      if bps == 2
        abit = qpsk_demod(rx_corr);
      end
      rx_bits = [rx_bits abit];
    end % c=2:Nc+1
    aphase_est_pilot_log = [aphase_est_pilot_log; aphase_est_pilot(2:Nc+1)];
  end 

  % Adjust nin to take care of sample clock offset

  nin = Nsamperframe;
  if timing_en && timing_valid
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

  % estimates of signal and noise power (see cohpsk.m for further explanation)
  % signal power is distance from axis on complex plane
  % we just measure noise power on imag axis, as it isn't affected by fading
  % using all symbols in frame worked better than just pilots
  
  sig_var = sum(abs(rx_np) .^ 2)/length(rx_np);
  sig_rms = sqrt(sig_var);
  
  sum_x = 0;
  sum_xx = 0;
  n = 0;
  for i=1:length(rx_np)
    s = rx_np(i);
    if abs(real(s)) > sig_rms 
      % select two constellation points on real axis
      sum_x  += imag(s);
      sum_xx += imag(s)*imag(s);
      n++;
    end
  end
   
  noise_var = 0;
  if n > 1
    noise_var = (n*sum_xx - sum_x*sum_x)/(n*(n-1));
  end

  % Total noise power is twice estimate of imaginary-axis noise.  This
  % effectively gives us the an estimate of Es/No
  
  states.noise_var = 2*noise_var; 
  states.sig_var = sig_var;

  % maintain mean amp estimate for LDPC decoder

  states.mean_amp = 0.9*states.mean_amp + 0.1*mean(rx_amp);
  
  states.rx_sym = rx_sym;
  states.rxbuf = rxbuf;
  states.nin = nin;
  states.timing_valid = timing_valid;
  states.timing_mx = timing_mx;
  states.timing_est = timing_est;
  states.sample_point = sample_point;
  states.delta_t = delta_t;
  states.foff_est_hz = foff_est_hz;
  states.coarse_foff_est_hz = coarse_foff_est_hz;
endfunction


% generate a test frame of ldpc encoded bits, used for raw and
% coded BER testing.  Includes UW and txt files

function [tx_bits payload_data_bits codeword] = create_ldpc_test_frame
  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 17; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;
  ldpc;
  gp_interleaver;
  
  % Set up LDPC code

  mod_order = 4; bps = 2; modulation = 'QPSK'; mapping = 'gray';

  init_cml('~/cml/'); % TODO: make this path sensible and portable
  load HRA_112_112.txt
  [code_param framesize rate] = ldpc_init_user(HRA_112_112, modulation, mod_order, mapping);
  assert(Nbitsperframe == (code_param.code_bits_per_frame + Nuwbits + Ntxtbits));

  rand('seed',1);
  payload_data_bits = round(rand(1,code_param.data_bits_per_frame));
  codeword = LdpcEncode(payload_data_bits, code_param.H_rows, code_param.P_matrix);
  Nsymbolsperframe = length(codeword)/bps;
  
  % need all these steps to get actual raw codeword bits at demod
  % note this will only work for single interleaver frame case,
  % but that's enough for initial quick tests
  
  tx_symbols = [];
  for s=1:Nsymbolsperframe
    tx_symbols = [tx_symbols qpsk_mod( codeword(2*(s-1)+1:2*s) )];
  end

  tx_symbols = gp_interleave(tx_symbols);

  codeword_raw = [];
  for s=1:Nsymbolsperframe
    codeword_raw = [codeword_raw qpsk_demod(tx_symbols(s))];
  end

  % insert UW and txt bits
  
  tx_bits = [tx_uw zeros(1,Ntxtbits) codeword_raw];
  assert(Nbitsperframe == length(tx_bits));

endfunction


% Save test bits frame to a text file in the form of a C array
% 
% usage:
%   ofdm_lib; test_bits_ofdm_file
%

function test_bits_ofdm_file
  [test_bits_ofdm payload_data_bits codeword] = create_ldpc_test_frame;
  printf("%d test bits\n", length(test_bits_ofdm));
  
  f=fopen("../src/test_bits_ofdm.h","wt");
  fprintf(f,"/* Generated by test_bits_ofdm_file() Octave function */\n\n");
  fprintf(f,"const int test_bits_ofdm[]={\n");
  for m=1:length(test_bits_ofdm)-1
    fprintf(f,"  %d,\n",test_bits_ofdm(m));
  endfor
  fprintf(f,"  %d\n};\n",test_bits_ofdm(end));

  fprintf(f,"\nconst int payload_data_bits[]={\n");
  for m=1:length(payload_data_bits)-1
    fprintf(f,"  %d,\n",payload_data_bits(m));
  endfor
  fprintf(f,"  %d\n};\n",payload_data_bits(end));

  fprintf(f,"\nconst int test_codeword[]={\n");
  for m=1:length(codeword)-1
    fprintf(f,"  %d,\n",codeword(m));
  endfor
  fprintf(f,"  %d\n};\n",codeword(end));

  fclose(f);

endfunction


% iterate state machine ------------------------------------

function states = sync_state_machine(states, rx_uw)
  ofdm_load_const;
  next_state = states.sync_state;
  states.sync_start = states.sync_end = 0;
  
  if strcmp(states.sync_state,'search') 

    if states.timing_valid

      % freq offset est has some bias, but this refinement step fixes bias

      st = M+Ncp + Nsamperframe + 1; en = st + 2*Nsamperframe;
      woff_est = 2*pi*states.foff_est_hz/Fs;
      [ct_est foff_est timing_valid timing_mx] = coarse_sync(states, states.rxbuf(st:en) .* exp(-j*woff_est*(st:en)), states.rate_fs_pilot_samples);
      if verbose
        printf("  coarse_foff: %4.1f refine: %4.1f combined: %4.1f\n", states.foff_est_hz, foff_est, states.foff_est_hz+foff_est);
      end
      states.foff_est_hz += foff_est;
      states.frame_count = 0;
      states.sync_counter = 0;
      states.sync_start = 1;
      next_state = 'trial';
    end
  end
        
  if strcmp(states.sync_state,'synced') || strcmp(states.sync_state,'trial')

    states.frame_count++;
    states.frame_count_interleaver++;
      
    % during trial sync we don't tolerate errors so much, once we have synced up
    % we are willing to wait out a fade
      
    if states.frame_count == 3
      next_state = 'synced';
    end
    if strcmp(states.sync_state,'synced')
      sync_counter_thresh = 12;
      uw_thresh = 2;
    else
      sync_counter_thresh = 2;
      uw_thresh = 1;
    end

    % freq offset est may be too far out, and has aliases every 1/Ts, so
    % we use a Unique Word to get a really solid indication of sync.

    states.uw_errors = sum(xor(tx_uw,rx_uw));
    if (states.uw_errors > uw_thresh)
      states.sync_counter++;
      if states.sync_counter == sync_counter_thresh
        next_state = 'search';
        states.sync_state_interleaver = 'search';
      end
    else
      states.sync_counter = 0;
    end
  end

  states.last_sync_state = states.sync_state;
  states.last_sync_state_interleaver = states.sync_state_interleaver;
  states.sync_state = next_state;
endfunction


% test function, kind of like a CRC for QPSK symbols, to compare two vectors

function acc = test_acc(v)
  sre = 0; sim = 0;
  for i=1:length(v)
    x = v(i);
    re = round(real(x)); im = round(imag(x));
    sre += re; sim += im;
    %printf("%d %10f %10f %10f %10f\n", i, re, im, sre, sim);
  end
  acc = sre + j*sim;
end
