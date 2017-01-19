% oqpsk.m
% David Rowe Jan 2017
%
% Unfiltered OQPSK modem implementation and simulations to test,
% derived from GMSK modem in gmsk.m

rand('state',1); 
randn('state',1);
graphics_toolkit ("gnuplot");
format

% init nodem states

function oqpsk_states = oqpsk_init(oqpsk_states, Rs)

  % general 

  verbose = oqpsk_states.verbose;
  oqpsk_states.Fs  = 4*Rs;
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

function [tx tx_symb] = oqpsk_mod(oqpsk_states, tx_bits)
  M = oqpsk_states.M;
  bps = oqpsk_states.bps;
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


#{ 

  Unfiltered OQPSK demodulator function, with (optional) phase and
  timing estimation.  Adapted from Fig 8 of [1].  See also gmsk.m and
  [2].

  Note demodulator returns phase corrected symbols sampled at ideal
  timing instant.  These symbols may have a m*pi/2 phase ambiguity due
  to properties of phase tracking loop.  teh aller is responsible for
  determining this ambiguity and recovering the actual bits.

  [1] GMSK demodulator in IEEE Trans on Comms, Muroyta et al, 1981,
  "GSM Modulation for Digital Radio Telephony".

  [2] GMSK Modem Simulation, http://www.rowetel.com/?p=3824

#}


function [rx_symb rx_int filt_log dco_log timing_adj Toff] = oqpsk_demod(oqpsk_states, rx)
  M = oqpsk_states.M;
  Rs = oqpsk_states.Rs;
  Fs = oqpsk_states.Fs;
  nsam = length(rx);
  nsym = floor(nsam/M);
  verbose = oqpsk_states.verbose;

  timing_angle_log = zeros(1,length(rx));
  rx_int = zeros(1,length(rx));
  dco_log = filt_log = zeros(1,nsam);

  % Unfiltered PSK - integrate energy in symbols M long in re and im arms

  rx_int = conv(rx,ones(1,M))/M;

  % phase and fine frequency tracking and correction ------------------------

  if oqpsk_states.phase_est
 
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
  end

  if  oqpsk_states.timing_est    
      % w is the ref sine wave at the timing clock frequency
      % tw is the length of the window used to estimate timing

      tw = 200*M;
      k = 1;
      xr_log = []; xi_log = [];
      w_log = [];
      timing_clock_phase = 0;
      timing_angle = pi;  % XXX
      timing_angle_log = zeros(1,nsam);
  end

  % Sample by sample processing loop for timing and phase est.  Note
  % this operates at sample rate Fs, unlike many PSK modems that
  % operate at the symbol rate Rs

  for i=1:nsam

    if oqpsk_states.timing_est

      % update sample timing estimate every tw samples, free wheel
      % rest of the time

      if mod(i,tw) == 0
        l = i - tw+1;
        xr = abs(real(rx_int(l:l+tw-1)));
        xi = abs(imag(rx_int(l:l+tw-1)));
        w = exp(j*(l:l+tw-1)*2*pi*Rs/Fs);
        X = xr * w';
        %timing_clock_phase = timing_angle = angle(X);
        timing_clock_phase = timing_angle = pi; % XXX
        k++;
        xr_log = [xr_log xr];
        xi_log = [xi_log xi];
        w_log = [w_log w];
      else
        timing_clock_phase += (2*pi)/M;
      end
      timing_angle_log(i) = timing_angle;
    end

    if oqpsk_states.phase_est

      % PLL per-sample processing

      rx_int(i) *= exp(-j*dco);
      ph_err = sign(real(rx_int(i))*imag(rx_int(i)))*cos(timing_clock_phase);
      lower = ph_err*G2 + lower;
      filt  = ph_err*G1 + lower;
      dco_log(i) = dco = pi; % XXX
      dco   = dco + filt;
      filt_log(i) = filt;

    end
  
  end

  % final adjustment of timing output to take into account slowly
  % moving estimates due to sample clock offset.  Unwrap ensures that
  % when timing angle jumps from -pi to pi we move to the next symbol
  % and frame sync isn't broken

  timing_adj = timing_angle_log*M/(2*pi); 
  timing_adj_uw = unwrap(timing_angle_log)*M/(2*pi); 
  % Toff = floor(2*M+timing_adj);
  Toff = floor(timing_adj_uw+0.5);

  % sample integrator output at correct timing instant
    
  k = 1;
  re_syms = im_syms = zeros(1,nsym);  
  rx_symb = [];
  for i=M:M:nsam
    if i-Toff(i)+M/2 <= nsam
      re_syms(k) = real(rx_int(i-Toff(i)));
      im_syms(k) = imag(rx_int(i-Toff(i)+M/2));
      %re_syms(k) = real(rx_int(i));
      %im_syms(k) = imag(rx_int(i));
      rx_symb = [rx_symb re_syms(k) + j*im_syms(k)];
      k++;
    end
  end
  
endfunction



% Test modem over a range Eb/No points in AWGN channel.  Doesn't resolve
% phase ambiguity or time offsets in received bit stream, so we can't test
% phase or timing offsets.  Does test basic lock of phase and timing loops.

function sim_out = oqpsk_test(sim_in)
  bitspertestframe =  sim_in.bitspertestframe;
  nbits =  sim_in.nbits;
  EbNodB = sim_in.EbNodB;
  verbose = sim_in.verbose;
  Rs = 4800;

  oqpsk_states.verbose = verbose;
  oqpsk_states.coherent_demod = sim_in.coherent_demod;
  oqpsk_states.phase_est     = sim_in.phase_est;
  oqpsk_states.timing_est    = sim_in.timing_est;
  oqpsk_states = oqpsk_init(oqpsk_states, Rs);
  M = oqpsk_states.M;
  Fs = oqpsk_states.Fs;
  Rs = oqpsk_states.Rs;
 
  tx_testframe = round(rand(1, bitspertestframe));
  ntestframes = floor(nbits/bitspertestframe);
  tx_bits = [];
  for i=1:ntestframes
    tx_bits = [tx_bits tx_testframe];
  end

  for ne = 1:length(EbNodB)
    aEbNodB = EbNodB(ne);
    EbNo = 10^(aEbNodB/10);
    variance = Fs/(Rs*EbNo*oqpsk_states.bps);

    [tx tx_symb] = oqpsk_mod(oqpsk_states, tx_bits);
    nsam = length(tx);
    
    noise = sqrt(variance/2)*(randn(1,nsam) + j*randn(1,nsam));
    st = 1+sim_in.timing_offset; en = length(tx);
    rx = tx(st:en)*exp(j*sim_in.phase_offset) + noise(st:en);

    [rx_symb rx_int filt_log dco_log timing_adj Toff] = oqpsk_demod(oqpsk_states, rx);
    
    % Determine ambiguities:
    %   a) could be m*pi/2 rotations of phase by phase est
    %   b) could be I and Q swapped by timing est
    %   c) time alignment of test frame

    nsymb = bitspertestframe/oqpsk_states.bps;
    nrx_symb = length(rx_symb);
    rx_bits = zeros(1, bitspertestframe);
    atx_symb = tx_symb(1:nsymb);

    % correlate with I and Q tx sequences at various offsets

    max_corr = real(atx_symb) * real(atx_symb)';
    for offset=2:nrx_symb-nsymb+1
      corr_ii(offset) = real(atx_symb) * real(rx_symb(offset:offset+nsymb-1))'/max_corr;
      corr_qq(offset) = imag(atx_symb) * imag(rx_symb(offset:offset+nsymb-1))'/max_corr;
      corr_iq(offset) = real(atx_symb) * imag(rx_symb(offset:offset+nsymb-1))'/max_corr;
      corr_qi(offset) = imag(atx_symb) * real(rx_symb(offset:offset+nsymb-1))'/max_corr;
      %printf("offset: %2d ii: % 5f qq: % 5f iq: % 5f qi: % 5f\n", 
      %offset, corr_ii(offset), corr_qq(offset), corr_iq(offset), corr_qi(offset));

      if abs(corr_ii(offset)) > 0.8
        % no IQ swap, or time offset
        arx_symb = real(rx_symb(offset:offset+nsymb-1)) + j*imag(rx_symb(offset:offset+nsymb-1));
      end
      if abs(corr_qi(offset)) > 0.8
        % IQ swap, I part in Q part of symbol before
        i_sign = sign(corr_iq(offset-1));
        q_sign = sign(corr_qi(offset));
        arx_symb = i_sign*imag(rx_symb(offset-1:offset+nsymb-2)) + j*q_sign*real(rx_symb(offset:offset+nsymb-1));
        for i=1:nsymb
          rx_bits(2*i-1:2*i) = qpsk_demod(arx_symb(i)*exp(-j*pi/4));
        end
        nerr = sum(xor(tx_testframe, rx_bits));
        printf("offset: %5d swap: %d i_sign: % 2.1f q_sign: % 2.1f nerr: %d\n", 
        offset, 1, i_sign, q_sign, nerr);
      end
    end


#{
    for i=1:nrx_symb-nsymb+1
      for k=0:3
        phase_amb = exp(j*k*pi/2);
        arx_symb = rx_symb(i:i+nsymb-1) .* phase_amb;
        for l=1:nsymb
          rx_bits(2*l-1:2*l) = qpsk_demod(arx_symb(l).*exp(-j*pi/4));
        end
        nerr = sum(xor(tx_testframe, rx_bits));
        if nerr == 0
          printf("i: %d k: %d nerr: %d\n", i, k, nerr);
        end
      end
    end
#}

    TERvec(ne) = 0;
    BERvec(ne) = 0;
#{
    Nerrs_min = nbits; Nbits_min = nbits; l = length(rx_bits);
    for i=1:1
      Nerrs = sum(xor(rx_bits(i:l), tx_bits(1:l-i+1)));
      if Nerrs < Nerrs_min
        Nerrs_min = Nerrs;
        Nbits_min = l;
      end
    end
 
    TERvec(ne) = Nerrs_min;
    BERvec(ne) = Nerrs_min/nbits;
    
    if verbose > 0
      printf("EbNo dB: %3.1f Nbits: %d Nerrs: %d BER: %4.3f BER Theory: %4.3f\n", 
      aEbNodB, nbits, Nerrs_min, BERvec(ne), 0.5*erfc(sqrt(EbNo)));
    end
#}

    figure(1); clf;
    subplot(211)
    stem(real(tx_symb))
    title('Tx symbols');
    subplot(212)
    stem(imag(tx_symb))

    figure(2); clf;
    f = fftshift(fft(rx));
    Tx = 20*log10(abs(f));
    plot((1:length(f))*Fs/length(f) - Fs/2, Tx)
    grid;
    title('OQPSK Demodulator Input Spectrum');

    figure(3); clf;
    nplot = min(16, nbits/oqpsk_states.bps);
    title('Rx Integrator');
    subplot(211)
    stem(real(rx_int(1:nplot*M)))
    axis([1 nplot*M -1 1])
    subplot(212)
    stem(imag(rx_int(1:nplot*M)))
    axis([1 nplot*M -1 1])

    figure(4); clf;
    subplot(211);
    plot(filt_log);
    title('PLL filter')
    subplot(212);
    plot(dco_log);
    title('PLL DCO phase');
    
    figure(5); clf;
    subplot(211)
    plot(timing_adj);
    title('Timing est');
    subplot(212)
    plot(Toff);
    title('Timing est unwrap');

    figure(6); clf;
    st = floor(0.5*nrx_symb);
    plot(rx_symb(st:nrx_symb), '+');
    title('Scatter Diagram');
    axis([-1.5 1.5 -1.5 1.5])

    figure(7); clf;
    subplot(211);
    stem(real(arx_symb));
    title('Rx symbols')
    subplot(212);
    stem(imag(arx_symb));
   
    figure(8); clf;
    subplot(211)
    stem(tx_testframe(1:min(20,length(rx_bits))))
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
  sim_in.coherent_demod   = 1;
  sim_in.phase_est        = 1;
  sim_in.timing_est       = 1;
  sim_in.bitspertestframe = 100;
  sim_in.nbits            = 1000;
  sim_in.EbNodB           = 100;
  sim_in.verbose          = 2;
  sim_in.phase_offset     = pi/2;
  sim_in.timing_offset    = 0;

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
  Rs = 4800;
  verbose = 1;
  aEbNodB = 10;
  phase_offset = 0;
  freq_offset  = 0;
  timing_offset = 1;
  sample_clock_offset_ppm = 0;
  nsym = Rs*10;

  oqpsk_states.verbose = verbose;
  oqpsk_states.coherent_demod = 1;
  oqpsk_states.phase_track    = 1;
  oqpsk_states = oqpsk_init(oqpsk_states, Rs);
  Fs = oqpsk_states.Fs;
  Rs = oqpsk_states.Rs;
  M  = oqpsk_states.M;

  % A frame consists of nsym random data bits.

  framesize = 480;
  nframes = floor(nsym/framesize);
  tx_frame = round(rand(1, framesize));
  tx_bits = [];
  for i=1:nframes
    tx_bits = [tx_bits tx_frame];
  end

  tx = oqpsk_mod(oqpsk_states, tx_bits);

  tx = resample(tx, 1E6, 1E6-sample_clock_offset_ppm);
  tx = [zeros(1,timing_offset) tx];
  nsam = length(tx);

  EbNo = 10^(aEbNodB/10);
  variance = Fs/(Rs*EbNo*oqpsk_states.bps);
  noise = sqrt(variance/2)*(randn(1,nsam) + j*randn(1,nsam));
  w  = (0:nsam-1)*2*pi*freq_offset/Fs + phase_offset;

  rx = sqrt(2)*tx.*exp(j*w) + noise;
  
  % optional dump to file

  if 0
    fc = 1500; gain = 10000;
    wc = 2*pi*fc/Fs;
    w1 = exp(j*wc*(1:nsam));
    rx1 = gain*real(rx .* w1);
    fout = fopen("rx_6dB.raw","wb");
    fwrite(fout, rx1, "short");
    fclose(fout);
  end

  % printf("ntx: %d nrx: %d ntx_bits: %d\n", length(tx), length(rx), length(tx_bits));

  [rx_bits rx_out rx_symb] = oqpsk_demod(oqpsk_states, rx);
  nframes_rx = length(rx_bits)/framesize;

  % printf("ntx: %d nrx: %d ntx_bits: %d nrx_bits: %d\n", length(tx), length(rx), length(tx_bits), length(rx_bits));

  % try different possible phase ambiguities, a UW would sort this out in the real world

  for i=1:4
    phase_amb = i*pi/2;
    rx_bits = [];
    for k=1:length(rx_symb)     
      rx_bits = [rx_bits qpsk_demod(rx_symb(k)*exp(j*(-pi/4+phase_amb)))];
    end

    [total_errors total_bits Nerrs_log Nerrs_all_log] = coarse_sync_ber(nframes_rx, tx_frame, rx_bits);

    if total_bits
      ber = total_errors/total_bits;
      printf("Eb/No: %3.1f f_off: %4.1f ph_off: %4.3f ph_amb: %4.3f Nframes: %d Nbits: %d Nerrs: %d BER: %f\n", 
             aEbNodB, freq_offset, phase_amb, phase_offset, nframes_rx, total_bits, total_errors, ber);

      figure(1); clf;
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

    end
  end

endfunction
    

% Choose one of these to run ------------------------------------------

run_oqpsk_single
%run_test_channel_impairments
%run_oqpsk_curves
%run_oqpsk_init

