% oqpsk.m
% David Rowe Jan 17
%
% Unfiltered OQPSK modem implementation and simulations to test,
% derived from GMSK

#{
  [ ] remove FM mod
  [ ] need OQPSK mod
  [ ] remove GMSK filters
  [ ] convert to QPSK from GMSK/BPSK 2T->T
  [ ] change names to oqpsk
#}

rand('state',1); 
randn('state',1);
graphics_toolkit ("gnuplot");
format

%
% Functions that implement the GMSK modem ------------------------------------------------------
%

function oqpsk_states = oqpsk_init(oqpsk_states, Rs)

  % general 

  verbose = oqpsk_states.verbose;
  oqpsk_states.Fs  = 48000;
  oqpsk_states.Rs  = Rs;
  oqpsk_states.bps = 2;              % two bit/symbol for QPSK    

  M = oqpsk_states.M = oqpsk_states.Fs/oqpsk_states.Rs;
  assert(floor(M) == M, "oversampling factor M must be an integer");
  assert(floor(M/2) == M/2, "(oversampling factor M)/2 must be an integer to offset QPSK");
endfunction


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
    if isscalar(symbol) == 0
        printf("only works with scalars\n");
        return;
    end
    bit0 = real(symbol*exp(j*pi/4)) < 0;
    bit1 = imag(symbol*exp(j*pi/4)) < 0;
    two_bits = [bit1 bit0];
endfunction


% Unfiltered OQPSK modulator 

function [tx tx_filt tx_symb] = oqpsk_mod(oqpsk_states, tx_bits)
  M = oqpsk_states.M;
  bps = oqpsk_states.bps
  nsym = length(tx_bits)/bps;
  nsam = nsym*M;
  verbose = oqpsk_states.verbose;

  % Map bits to Gray coded QPSK symbols

  tx_symb = zeros(1,nsym);
  for i=1:nsym
    tx_symb(i) = qpsk_mod(tx_bits(2*i-1:2*i))*exp(j*pi/4);
  end
  
  % Oversample by M (sample and hold) to create unfiltered QPSK

  tx = zeros(1, nsam);
  for i=1:nsym
    tx((i-1)*M+1:(i*M)) = tx_symb(i);
  end
  
  % delay Q arm by half of a symbol to make OQPSK

  tx = [real(tx) zeros(1,M/2)] + j*[zeros(1,M/2) imag(tx)];
endfunction


function [rx_bits rx_int rx_symb] = oqpsk_demod(oqpsk_states, rx)
  M = oqpsk_states.M;
  Rs = oqpsk_states.Rs;
  Fs = oqpsk_states.Fs;
  nsam = length(rx);
  nsym = floor(nsam/M);

  timing_angle_log = zeros(1,length(rx));
  rx_int = zeros(1,length(rx));

  % Unfiltered PSK - integrate energy in symbols M long in re and im arms

  rx_int = conv(rx,ones(1,M))/M;

  % phase and fine frequency tracking and correction ------------------------

  if oqpsk_states.phase_track
 
      % DCO design from "Introduction To Phase-Lock Loop System Modeling", Wen Li
      % http://www.ece.ualberta.ca/~ee401/parts/data/PLLIntro.pdf

      eta = 0.707;
      wn = 2*pi*10*(Rs/4800);  % (Rs/4800) -> found reducing the BW benefical with falling Rs
      Ts = 1/Fs;
      g1 = 1 - exp(-2*eta*wn*Ts);
      g2 = 1 + exp(-2*eta*wn*Ts) - 2*exp(-eta*wn*Ts)*cos(wn*Ts*sqrt(1-eta*eta));
      Gpd = 2/pi;
      Gvco = 1;
      G1 = g1/(Gpd*Gvco);  G2 = g2/(Gpd*Gvco);
      %printf("g1: %e g2: %e G1: %e G2: %e\n", g1, g2, G1, G2);

      filt_prev = dco = lower = ph_err_filt = ph_err = 0;
      dco_log = filt_log = zeros(1,nsam);

      % w is the ref sine wave at the timing clock frequency
      % tw is the length of the window used to estimate timing

      k = 1;
      tw = 200*M;
      xr_log = []; xi_log = [];
      w_log = [];
      timing_clock_phase = 0;
      timing_angle = 0;
      timing_angle_log = zeros(1,nsam);

      for i=1:nsam

        % update sample timing estimate every tw samples

        if mod(i,tw) == 0
          l = i - tw+1;
          xr = abs(real(rx_int(l:l+tw-1)));
          xi = abs(imag(rx_int(l:l+tw-1)));
          w = exp(j*(l:l+tw-1)*2*pi*(Rs/2)/Fs);
          X = xr * w';
          timing_clock_phase = timing_angle = angle(X);
          k++;
          xr_log = [xr_log xr];
          xi_log = [xi_log xi];
          w_log = [w_log w];
        else
          timing_clock_phase += (2*pi)/(2*M);
        end
        timing_angle_log(i) = timing_angle;

        rx_int(i) *= exp(-j*dco);
        ph_err = sign(real(rx_int(i))*imag(rx_int(i)))*cos(timing_clock_phase);
        lower = ph_err*G2 + lower;
        filt  = ph_err*G1 + lower;
        dco   = dco + filt;
        filt_log(i) = filt;
        dco_log(i) = dco;
      end
      
      figure(10);
      clf
      subplot(211);
      plot(filt_log);
      title('PLL filter')
      subplot(212);
      plot(dco_log/pi);
      title('PLL DCO phase');
      %axis([1 nsam -0.5 0.5])
  end

  % sample integrator output at correct timing instant
    
  timing_adj = timing_angle_log*2*M/(2*pi);
  timing_adj_uw = unwrap(timing_angle_log)*2*M/(2*pi);
  % Toff = floor(2*M+timing_adj);
  Toff = floor(timing_adj_uw+0.5);

  k = 1;
  re_syms = im_syms = zeros(1,nsym);  
  rx_bits = []; rx_symb = [];
  for i=M:M:nsam
    re_syms(k) = real(rx_int(i-Toff(i)));
    im_syms(k) = imag(rx_int(i-Toff(i)+M/2));
    %re_syms(k) = real(rx_int(i));
    %im_syms(k) = imag(rx_int(i));
    rx_symb = [rx_symb re_syms(k) + j*im_syms(k)];
    rx_bits = [rx_bits qpsk_demod((re_syms(k) + j*im_syms(k))*exp(-j*pi/4))];
    k++;
  end
  
  
  figure(11);
  subplot(211)
  plot(timing_adj);
  title('Timing est');
  subplot(212)
  plot(Toff);
  title('Timing est unwrap');

endfunction



%  Functions for Testing the OQPSK modem --------------------------------------------------------
%

function sim_out = oqpsk_test(sim_in)
  nbits =  sim_in.nbits;
  EbNodB = sim_in.EbNodB;
  verbose = sim_in.verbose;
  Rs = 4800;

  oqpsk_states.verbose = verbose;
  oqpsk_states.coherent_demod = sim_in.coherent_demod;
  oqpsk_states.phase_track    = sim_in.phase_track;
  oqpsk_states = oqpsk_init(oqpsk_states, Rs);
  M = oqpsk_states.M;
  Fs = oqpsk_states.Fs;
  Rs = oqpsk_states.Rs;
 
  for ne = 1:length(EbNodB)
    aEbNodB = EbNodB(ne);
    EbNo = 10^(aEbNodB/10);
    variance = Fs/(Rs*EbNo*oqpsk_states.bps);

    tx_bits = round(rand(1, nbits));
    %tx_bits = zeros(1,nbits);
    [tx tx_filt tx_symbols] = oqpsk_mod(oqpsk_states, tx_bits);
    nsam = length(tx);
    
    noise = sqrt(variance/2)*(randn(1,nsam) + j*randn(1,nsam));
    rx    = tx + noise;

    [rx_bits rx_out rx_symb] = oqpsk_demod(oqpsk_states, rx(1:length(rx)));

    % search for frame location over a range

    Nerrs_min = nbits; Nbits_min = nbits; l = length(rx_bits);
    for i=1:1
      Nerrs = sum(xor(rx_bits(i:l), tx_bits(1:l-i+1)));
      if Nerrs < Nerrs_min
        Nerrs_min = Nerrs;
        Nbits_min = l;
      end
    end
 
    TERvec(ne) = Nerrs_min;
    BERvec(ne) = Nerrs_min/nbits
    
    if verbose > 0
      printf("EbNo dB: %3.1f Nbits: %d Nerrs: %d BER: %4.3f BER Theory: %4.3f\n", 
      aEbNodB, nbits, Nerrs_min, BERvec(ne), 0.5*erfc(sqrt(EbNo)));
    end

    figure(1); clf;
    subplot(211)
    plot(real(tx))
    subplot(212)
    plot(imag(tx))
    title('OQPSK tx sequence');

    figure(2); clf;
    f = fftshift(fft(rx));
    Tx = 20*log10(abs(f));
    plot(Tx)
    grid;
    title('OQPSK Demodulator Input Spectrum');

    figure(3); clf;
    nplot = min(16, nbits/oqpsk_states.bps);
    subplot(211)
    plot(real(rx_out(1:nplot*M))/(M))
    title('Integrator');
    axis([1 nplot*M -1 1])
    subplot(212)
    plot(imag(rx_out(1:nplot*M)/(M)))
    axis([1 nplot*M -1 1])
    title('Rx integrator');

    figure(3); clf;
    plot(rx_symb, '+');
    title('Scatter Diagram');

    figure(4); clf;
    subplot(211)
    stem(tx_bits(1:min(20,nbits)))
    title('Tx Bits')
    subplot(212)
    stem(rx_bits(1:min(20,length(rx_bits))))
    title('Rx Bits')
  end

  sim_out.TERvec = TERvec;
  sim_out.BERvec = BERvec;
  sim_out.Rs = oqpsk_states.Rs;
endfunction


function run_oqpsk_single
  sim_in.coherent_demod = 1;
  sim_in.phase_track    = 0;
  sim_in.nbits          = 10000;
  sim_in.EbNodB         = 4;
  sim_in.verbose        = 2;

  sim_out = oqpsk_test(sim_in);
endfunction


% Generate a bunch of BER versus Eb/No curves for various demods

function run_oqpsk_curves
  sim_in.coherent_demod = 1;
  sim_in.nsym = 48000;
  sim_in.EbNodB = 2:10;
  sim_in.verbose = 1;

  oqpsk_coh = oqpsk_test(sim_in);

  sim_in.coherent_demod = 0;
  oqpsk_noncoh = oqpsk_test(sim_in);

  Rs = oqpsk_coh.Rs;
  EbNo  = 10 .^ (sim_in.EbNodB/10);
  alpha = 0.75; % guess for BT=0.5 OQPSK
  oqpsk_theory.BERvec = 0.5*erfc(sqrt(alpha*EbNo));

  % BER v Eb/No curves

  figure;
  clf;
  semilogy(sim_in.EbNodB, oqpsk_theory.BERvec,'r;OQPSK theory;')
  hold on;
  semilogy(sim_in.EbNodB, oqpsk_coh.BERvec,'g;OQPSK sim coherent;')
  semilogy(sim_in.EbNodB, oqpsk_noncoh.BERvec,'b;OQPSK sim non-coherent;')
  hold off;
  grid("minor");
  axis([min(sim_in.EbNodB) max(sim_in.EbNodB) 1E-4 1])
  legend("boxoff");
  xlabel("Eb/No (dB)");
  ylabel("Bit Error Rate (BER)")

  % BER v C/No (1 Hz noise BW and Eb=C/Rs=1/Rs)
  % Eb/No = (C/Rs)/(1/(N/B))
  % C/N   = (Eb/No)*(Rs/B)

  RsOnB_dB = 10*log10(Rs/1);
  figure;
  clf;
  semilogy(sim_in.EbNodB+RsOnB_dB, oqpsk_theory.BERvec,'r;OQPSK theory;')
  hold on;
  semilogy(sim_in.EbNodB+RsOnB_dB, oqpsk_coh.BERvec,'g;OQPSK sim coherent;')
  semilogy(sim_in.EbNodB+RsOnB_dB, oqpsk_noncoh.BERvec,'b;OQPSK sim non-coherent;')
  hold off;
  grid("minor");
  axis([min(sim_in.EbNodB+RsOnB_dB) max(sim_in.EbNodB+RsOnB_dB) 1E-4 1])
  legend("boxoff");
  xlabel("C/No for Rs=4800 bit/s and 1 Hz noise bandwidth (dB)");
  ylabel("Bit Error Rate (BER)")

endfunction



% attempt to perform "coarse sync" sync with the received frames, we
% check each frame for the best coarse sync position.  Brute force
% approach, that would be changed for a real demod which has some
% sort of unique word.  Start looking for valid frames 1 frame
% after start of pre-amble to give PLL time to lock

function [total_errors total_bits Nerrs_log Nerrs_all_log errors_log] = coarse_sync_ber(nframes_rx, tx_frame, rx_bits)

  Nerrs_log = zeros(1, nframes_rx);
  Nerrs_all_log = zeros(1, nframes_rx);
  total_errors = 0;
  total_bits   = 0;
  framesize = length(tx_frame);
  errors_log = [];

  for f=2:nframes_rx-1
    Nerrs_min = framesize;
    for i=1:framesize;
      st = (f-1)*framesize+i; en = st+framesize-1;
      errors = xor(rx_bits(st:en), tx_frame); 
      Nerrs = sum(errors);
      if Nerrs < Nerrs_min
        Nerrs_min = Nerrs;
        errors_min = errors;
      end
    end
    Nerrs_all_log(f) = Nerrs_min;
    if Nerrs_min/framesize < 0.1
      errors_log = [errors_log errors_min];
      Nerrs_log(f) = Nerrs_min;
      total_errors += Nerrs_min;
      total_bits   += framesize;
    end
  end
endfunction


function plot_spectrum(oqpsk_states, rx, preamble_location, title_str)
  Fs = oqpsk_states.Fs;
  st = preamble_location + oqpsk_states.npreamble*oqpsk_states.M;
  sig = rx(st:st+Fs*0.5);
  h = hanning(length(sig))';
  Rx=20*log10(abs(fftshift(fft(sig .* h, Fs))));
  figure;
  plot(-Fs/2:Fs/2-1,Rx);
  grid("minor");
  xlabel('Hz');
  ylabel('dB');
  topy = ceil(max(Rx)/10)*10;
  axis([-4000 4000 topy-50 topy+10])
  title(title_str);
endfunction

% Give the demod a hard time: frequency, phase, time offsets, sample clock difference
   
function run_test_channel_impairments
  Rs = 1200;
  verbose = 1;
  aEbNodB = 6;
  phase_offset = pi/2;
  freq_offset  = -104;
  timing_offset = 100E3;
  sample_clock_offset_ppm = -500;
  interferer_freq = -1500;
  interferer_amp  = 0;
  nsym = 4800*2;
  npreamble = 480;

  oqpsk_states.npreamble = npreamble;
  oqpsk_states.verbose = verbose;
  oqpsk_states.coherent_demod = 1;
  oqpsk_states.phase_track    = 1;
  oqpsk_states = oqpsk_init(oqpsk_states, Rs);
  Fs = oqpsk_states.Fs;
  Rs = oqpsk_states.Rs;
  M  = oqpsk_states.M;

  % A frame consists of nsym random data bits.  Some experimentation
  % has shown they must be random-ish data (not say 11001100...) for
  % timing estimator to work.  However initial freq offset estimation
  % is a lot easier with a 01010 type sequence, so we construct a 
  % frame with a pre-amble followed by frames of random data.

  framesize = 480;
  nframes = floor(nsym/framesize);
  tx_frame = round(rand(1, framesize));
  tx_bits = zeros(1,npreamble);
  tx_bits(1:2:npreamble) = 1;
  for i=1:nframes
    tx_bits = [tx_bits tx_frame];
  end

  [tx tx_filt tx_symbols] = oqpsk_mod(oqpsk_states, tx_bits);

  tx = resample(tx, 1E6, 1E6-sample_clock_offset_ppm);
  tx = [zeros(1,timing_offset) tx];
  nsam = length(tx);

  if verbose > 1
    figure;
    subplot(211)
    st = timing_offset; en = st+M*10;
    plot(real(tx(st:en)))
    title('Real part of tx');
    subplot(212)
    plot(imag(tx(st:en)))
    title('Imag part of tx');
  end

  EbNo = 10^(aEbNodB/10);
  variance = Fs/(Rs*EbNo);
  noise = sqrt(variance/2)*(randn(1,nsam) + j*randn(1,nsam));
  w  = (0:nsam-1)*2*pi*freq_offset/Fs + phase_offset;
  interferer = interferer_amp*exp(j*interferer_freq*(2*pi/Fs)*(0:nsam-1));

  rx = sqrt(2)*tx.*exp(j*w) + noise + interferer;
  
  % optional dump to file

  if 1
    fc = 1500; gain = 10000;
    wc = 2*pi*fc/Fs;
    w1 = exp(j*wc*(1:nsam));
    rx1 = gain*real(rx .* w1);
    fout = fopen("rx_6dB.raw","wb");
    fwrite(fout, rx1, "short");
    fclose(fout);
  end

  rx = rx1 .* conj(w1);

  [preamble_location freq_offset_est] = find_preamble(oqpsk_states, M, npreamble, rx);
  w_est  = (0:nsam-1)*2*pi*freq_offset_est/Fs;
  rx = rx.*exp(-j*w_est);

  plot_spectrum(oqpsk_states, rx, preamble_location, "OQPSK rx just after preamble");

  % printf("ntx: %d nrx: %d ntx_bits: %d\n", length(tx), length(rx), length(tx_bits));

  [rx_bits rx_out rx_filt] = oqpsk_demod(oqpsk_states, rx(preamble_location+framesize:nsam));
  nframes_rx = length(rx_bits)/framesize;

  % printf("ntx: %d nrx: %d ntx_bits: %d nrx_bits: %d\n", length(tx), length(rx), length(tx_bits), length(rx_bits));

  [total_errors total_bits Nerrs_log Nerrs_all_log] = coarse_sync_ber(nframes_rx, tx_frame, rx_bits);

  ber = total_errors/total_bits;

  printf("Eb/No: %3.1f f_off: %4.1f ph_off: %4.3f Nframes: %d Nbits: %d Nerrs: %d BER: %f\n", 
         aEbNodB, freq_offset, phase_offset, nframes_rx, total_bits, total_errors, ber);

  figure;
  clf
  subplot(211)
  plot(Nerrs_log,'r;errors/frame counted for BER;');
  hold on;
  plot(Nerrs_all_log,'g;all errors/frame;');
  hold off;
  legend("boxoff");
  title('Bit Errors')
  subplot(212)
  stem(real(cumsum(Nerrs_log)))
  title('Cumulative Bit Errors')

endfunction
    

% Choose one of these to run ------------------------------------------

run_oqpsk_single
%run_oqpsk_curves
%run_oqpsk_init

