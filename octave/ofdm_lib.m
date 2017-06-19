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


% Correlates the OFDM pilot symbol samples with a window of received
% samples to determine the most likely timing offset.  Combines two
% frames pilots so we need at least Nsamperframe+M+Ncp samples in rx.
% Also determines frequency offset at maximimum correlation.  Can be
% used for acquisition (coarse timing a freq offset), and fine timing

#{
  TODO: 
    [ ] attempt to speed up sync
        + tis rather stateless, this current demod, which is nice
        [ ] 0.5Hz grid and measure BER
            + so run demod a bunch of times at different offsets
            + Hmm cld also try +/- 20Hz multiples as it's aliased?
            + might to use 
            + need metric for sync, could be callback.
        [ ] or refine freq offset using pilots
        [ ] different error measure that 10% maybe soft dec
            + 10% very high BER
    [ ] simpler CPU/DFT for freq offset estimation
        + more suitable for real time implementation
#}

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

  states.foff_est_gain = 0.01;
  states.foff_est_hz = 0;
  states.sample_point = states.timing_est = 1;
  states.nin = states.Nsamperframe;

  % generate OFDM pilot symbol, used for timing and freq offset est

  rate_fs_pilot_samples = states.pilots * W/states.M;
  states.rate_fs_pilot_samples = [rate_fs_pilot_samples(states.M-states.Ncp+1:states.M) rate_fs_pilot_samples];
  
  % LDPC code is optionally enabled

  states.rate = 1.0;
  states.ldpc_en = 0;

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

  delta_t = 0;
  if timing_en
    % update timing at start of every frame

    st = M+Ncp + Nsamperframe + 1 - floor(ftwindow_width/2) + (timing_est-1);
    en = st + Nsamperframe-1 + M+Ncp + ftwindow_width-1;
          
    ft_est = coarse_sync(states, rxbuf(st:en) .* exp(-j*woff_est*(st:en)), rate_fs_pilot_samples);
    timing_est = timing_est + ft_est - ceil(ftwindow_width/2);

    if verbose > 1
      printf("  ft_est: %2d timing_est: %2d sample_point: %2d\n", ft_est, timing_est, sample_point);
    end

    % Black magic to keep sample_point inside cyclic prefix.  Or something like that.

    delta_t = ft_est - ceil(ftwindow_width/2);
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

  states.rx_sym = rx_sym;
  states.rxbuf = rxbuf;
  states.nin = nin;
  states.timing_est = timing_est;
  states.sample_point = sample_point;
  states.delta_t = delta_t;
  states.foff_est_hz = foff_est_hz;
endfunction


% Save test bits frame to a text file in the form of a C array
% 
% usage:
%   ofdm_lib; test_bits_ofdm_file
%

function test_bits_ofdm_file

  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 16; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;

  rand('seed',1);
  test_bits_ofdm = round(rand(1,Nbitsperframe));

  f=fopen("../src/test_bits_ofdm.h","wt");
  fprintf(f,"/* Generated by test_bits_ofdm_file() Octave function */\n\n");
  fprintf(f,"const int test_bits_ofdm[]={\n");
  for m=1:length(test_bits_ofdm)-1
    fprintf(f,"  %d,\n",test_bits_ofdm(m));
  endfor
  fprintf(f,"  %d\n};\n",test_bits_ofdm(length(test_bits_ofdm)));
  fclose(f);

endfunction
