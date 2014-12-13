% fsk.m
% David Rowe Nov 2014

% Simulation to test FSK demod
%
% TODO
%   [X] Code up mod/non-coh demod/AWGN channel simulation
%   [X] Eb/No verses BER curves
%   [X] test analog FM with pre/de-emphahsis 
%       + this will introduce delay, use fir filter, group delay
%   [X] channel simulation of HT/FM radio
%       + filtering, varying modulation index
%   [ ] GMSK
%   [X] refactor to plot analog FM demod curves
%   [X] SSB curves
%   [-] different modn index beta curves, 
%       + not really important for now
%   [ ] plot to illustrate harmonic dist, 
%   [X] integration with AFSK, AFSK-FM, v AFSK-SSB (ie FSK)
%   [-] fine timing offset for pre/de filters?
%       + do we need interpolation as well?
%       + might leave this as pre/de not significant now
%   [X] C/No curves?
%   [ ] spectrum plots

rand('state',1); 
randn('state',1);
graphics_toolkit ("gnuplot");

function sim_out = fsk_ber_test(sim_in)
  Fs        = 96000;
  fmark     = sim_in.fmark;
  fspace    = sim_in.fspace;
  Rs        = sim_in.Rs;
  Ts        = Fs/Rs;
  emphasis  = 50E-6;
  verbose   = sim_in.verbose;

  nsym      = sim_in.nsym;
  nsam      = nsym*Ts;
  EbNodB    = sim_in.EbNodB;

  fm        = sim_in.fm;

  if fm
    fm_states.pre_emp = 0;
    fm_states.de_emp  = 0;
    fm_states.Ts      = Ts;
    fm_states = analog_fm_init(fm_states);
  end

  % simulate over a range of Eb/No values

  for ne = 1:length(EbNodB)
    Nerrs = Terrs = Tbits = 0;

    % randn() generates noise across the entire Fs bandwidth, we want to scale
    % the noise power (i.e. the variance) so we get the Eb/No we want in the
    % bandwidth of our FSK signal, which we assume is Rs.

    aEbNodB = EbNodB(ne);
    EbNo = 10^(aEbNodB/10);
    variance = Fs/(Rs*EbNo);

    % Modulator -------------------------------

    tx_bits = round(rand(1, nsym));
    tx = zeros(1,nsam);
    tx_phase = 0;

    for i=1:nsym
      for k=1:Ts
        if tx_bits(i) == 1
          tx_phase += 2*pi*fmark/Fs;
        else
          tx_phase += 2*pi*fspace/Fs;
        end
        tx_phase = tx_phase - floor(tx_phase/(2*pi))*2*pi;
        tx((i-1)*Ts+k) = exp(j*tx_phase);
      end
    end

    % Optional AFSK over FM modulator

    if sim_in.fm
      % FM mod takes real input; +/- 1 for correct deviation
      tx = analog_fm_mod(fm_states, real(tx));
    end
  
    % Channel ---------------------------------

    % We use complex (single sided) channel simulation, as it's convenient
    % for the FM simulation.

    noise = sqrt(variance/2)*(randn(1,nsam) + j*randn(1,nsam));
    rx    = tx + noise;
    if verbose > 1
      printf("EbNo: %f Eb: %f var No: %f EbNo (meas): %f\n", 
      EbNo, var(tx)*Ts/Fs, var(noise)/Fs, (var(tx)*Ts/Fs)/(var(noise)/Fs));
    end

    % Optional AFSK over FM demodulator

    if sim_in.fm
      % scaling factor for convenience to match pure FSK
      rx_bb = 2*analog_fm_demod(fm_states, rx);
    else
      rx_bb = rx;
    end

    % Demodulator -----------------------------

    % non-coherent FSK demod

    mark_dc  = rx_bb .* exp(-j*(0:nsam-1)*2*pi*fmark/Fs);
    space_dc = rx_bb .* exp(-j*(0:nsam-1)*2*pi*fspace/Fs);

    rx_bits = zeros(1, nsym);
    for i=1:nsym
      st = (i-1)*Ts+1;
      en = st+Ts-1;
      mark_int(i)  = sum(mark_dc(st:en));
      space_int(i) = sum(space_dc(st:en));
      rx_bits(i) = abs(mark_int(i)) > abs(space_int(i));
    end
  
    if fm
      d = fm_states.nsym_delay;
      error_positions = xor(rx_bits(1+d:nsym), tx_bits(1:(nsym-d)));
    else
      error_positions = xor(rx_bits, tx_bits);
    end
    Nerrs = sum(error_positions);
    Terrs += Nerrs;
    Tbits += length(error_positions);

    TERvec(ne) = Terrs;
    BERvec(ne) = Terrs/Tbits;

    if verbose > 1
      figure(2)
      clf
      Rx = 10*log10(abs(fft(rx)));
      plot(Rx(1:Fs/2));
      axis([1 Fs/2 0 50]);

      figure(3)
      clf;
      subplot(211)
      plot(real(rx_bb(1:Ts*20)))
      subplot(212)
      Rx_bb = 10*log10(abs(fft(rx_bb)));
      plot(Rx_bb(1:3000));
      axis([1 3000 0 50]);

      figure(4);
      subplot(211)
      stem(abs(mark_int(1:100)));
      subplot(212)
      stem(abs(space_int(1:100)));   
      
      figure(5)
      clf
      plot(error_positions);
    end

    if verbose
      printf("EbNo (db): %3.2f Terrs: %d BER: %3.2f \n", aEbNodB, Terrs, Terrs/Tbits);
    end
  end

  sim_out.TERvec = TERvec;
  sim_out.BERvec = BERvec;
endfunction


function fm_states = analog_fm_init(fm_states)

  % FM modulator constants

  fm_states.Fs = Fs = 96000; FsOn2 = Fs/2;  
  fm_states.fm_max = fm_max = 3000;          % max modulation freq
  fm_states.fc = 24E3;                       % carrier frequency
  fm_states.fd = fd = 5E3;                   % (max) deviation
  fm_states.m = fd/fm_max;                   % modulation index
  fm_states.Bfm = Bfm = 2*(fd+fm_max);       % Carson's rule for FM signal bandwidth
  fm_states.tc = tc = 50E-6;
  fm_states.prede = [1 -(1 - 1/(tc*Fs))];    % pre/de emp filter coeffs

  % Select length of filter to be an integer number of symbols to
  % assist with "fine" timing offset estimation.  Set Ts to 1 for
  % analog modulation.

  Ts = fm_states.Ts;
  desired_ncoeffs = 200;
  ncoeffs = floor(desired_ncoeffs/Ts+1)*Ts;

  % "coarse" timing offset is half filter length, we have two filters.
  % This is the delay the two filters introduce, so we need to adjust
  % for this when comparing tx to trx bits for BER calcs.

  fm_states.nsym_delay = ncoeffs/Ts;

  % input filter gets rid of excess noise before demodulator, as too much
  % noise causes atan2() to jump around, e.g. -pi to pi.  However this
  % filter can cause harmonic distortion at very high SNRs, as it knocks out
  % some of the FM signal spectra.  This filter isn't really required for high
  % SNRs > 20dB.

  fc = (Bfm/2)/(FsOn2);
  fm_states.bin  = firls(ncoeffs,[0 fc*(1-0.05) fc*(1+0.05) 1],[1 1 0.01 0.01]);

  % demoduator output filter to limit us to fm_max (e.g. 3kHz)

  fc = (fm_max)/(FsOn2);
  fm_states.bout = firls(ncoeffs,[0 0.95*fc 1.05*fc 1], [1 1 0.01 0.01]);

endfunction


function tx = analog_fm_mod(fm_states, mod)
  Fs = fm_states.Fs;
  fc = fm_states.fc; wc = 2*pi*fc/Fs;
  fd = fm_states.fd; wd = 2*pi*fd/Fs;
  nsam = length(mod);

  if fm_states.pre_emp
    mod = filter(fm_states.prede,1,mod);
    mod = mod/max(mod);           % AGC to set deviation
  end

  tx_phase = 0;
  tx = zeros(1,nsam);

  for i=0:nsam-1
    w = wc + wd*mod(i+1);
    tx_phase = tx_phase + w;
    tx_phase = tx_phase - floor(tx_phase/(2*pi))*2*pi;
    tx(i+1) = exp(j*tx_phase);
  end       
endfunction


function [rx_out rx_bb] = analog_fm_demod(fm_states, rx)
  Fs = fm_states.Fs;
  fc = fm_states.fc; wc = 2*pi*fc/Fs;
  fd = fm_states.fd; wd = 2*pi*fd/Fs;
  nsam = length(rx);
  t = 0:(nsam-1);

  rx_bb = rx .* exp(-j*wc*t);      % down to complex baseband
  rx_bb = filter(fm_states.bin,1,rx_bb);
  rx_bb_diff = [ 1 rx_bb(2:nsam) .* conj(rx_bb(1:nsam-1))];
  rx_out = (1/wd)*atan2(imag(rx_bb_diff),real(rx_bb_diff));
  rx_out = filter(fm_states.bout,1,rx_out);
  if fm_states.de_emp
    rx_out = filter(1,fm_states.prede,rx_out);
  end
endfunction


function sim_out = analog_fm_test(sim_in)
  nsam      = sim_in.nsam;
  CNdB      = sim_in.CNdB;
  verbose   = sim_in.verbose;

  fm_states.pre_emp = pre_emp = sim_in.pre_emp;
  fm_states.de_emp  = de_emp = sim_in.de_emp;
  fm_states.Ts = 1;
  fm_states = analog_fm_init(fm_states);
  sim_out.Bfm = fm_states.Bfm;

  Fs = fm_states.Fs;
  Bfm = fm_states.Bfm;
  m = fm_states.m; tc = fm_states.tc; fm_max = fm_states.fm_max;
  t = 0:(nsam-1);

  fm = 1000; wm = 2*pi*fm/fm_states.Fs;
  
  % start simulation

  for ne = 1:length(CNdB)

    % work out the variance we need to obtain our C/N in the bandwidth
    % of the FM demod.  The gaussian generator randn() generates noise
    % with a bandwidth of Fs

    aCNdB = CNdB(ne);
    CN = 10^(aCNdB/10);
    variance = Fs/(CN*Bfm);
     
    % FM Modulator -------------------------------

    mod = sin(wm*t);
    tx = analog_fm_mod(fm_states, mod);

    % Channel ---------------------------------

    noise = sqrt(variance/2)*(randn(1,nsam) + j*randn(1,nsam));
    rx = tx + noise;

    % FM Demodulator

    [rx_out rx_bb] = analog_fm_demod(fm_states, rx);

    % notch out test tone

    w = 2*pi*fm/Fs; beta = 0.99;
    rx_notch = filter([1 -2*cos(w) 1],[1 -2*beta*cos(w) beta*beta], rx_out);

    % measure power with and without test tone to determine S+N and N

    settle = 1000;             % filter settling time, to avoid transients
    nsettle = nsam - settle;

    sinad = (rx_out(settle:nsam) * rx_out(settle:nsam)')/nsettle;
    nad = (rx_notch(settle:nsam) * rx_notch(settle:nsam)')/nsettle;

    snr = (sinad-nad)/nad;
    sim_out.snrdB(ne) = 10*log10(snr);
   
    % Theory from FMTutorial.pdf, Lawrence Der, Silicon labs paper

    snr_theory_dB = aCNdB + 10*log10(3*m*m*(m+1));
    fx = 1/(2*pi*tc); W = fm_max;
    I = (W/fx)^3/(3*((W/fx) - atan(W/fx)));
    I_dB = 10*log10(I);

    sim_out.snr_theorydB(ne) = snr_theory_dB;
    sim_out.snr_theory_pre_dedB(ne) = snr_theory_dB + I_dB;
   
    if verbose > 1
      printf("modn index: %2.1f Bfm: %.0f Hz\n", m, Bfm);
    end

    if verbose > 0
      printf("C/N: %4.1f SNR: %4.1f dB THEORY: %4.1f dB or with pre/de: %4.1f dB\n", 
      aCNdB, 10*log10(snr), snr_theory_dB, snr_theory_dB+I_dB);
    end

    if verbose > 1
      figure(1)
      subplot(211)
      plot(20*log10(abs(fft(rx))))
      title('FM Modulator Output Spectrum');
      axis([1 length(tx) 0 100]);
      subplot(212)
      Rx_bb = 20*log10(abs(fft(rx_bb)));
      plot(Rx_bb)
      axis([1 length(tx) 0 100]);
      title('FM Demodulator (baseband) Input Spectrum');

      figure(2)
      subplot(211)
      plot(rx_out(settle:nsam))
      axis([1 4000 -1 1]) 
      subplot(212)
      Rx = 20*log10(abs(fft(rx_out(settle:nsam))));
      plot(Rx(1:10000))
      axis([1 10000 0 100]);
   end

  end

endfunction

more off;

function run_fsk_curves
  sim_in.fmark     = 1200;
  sim_in.fspace    = 2200;
  sim_in.Rs        = 1200;
  sim_in.nsym      = 12000;
  sim_in.EbNodB    = 0:2:20;
  sim_in.fm        = 0;
  sim_in.verbose   = 1;

  EbNo  = 10 .^ (sim_in.EbNodB/10);
  fsk_theory.BERvec = 0.5*exp(-EbNo/2); % non-coherent BFSK demod
  fsk_sim = fsk_ber_test(sim_in);

  sim_in.fm  = 1;
  fsk_fm_sim = fsk_ber_test(sim_in);

  % BER v Eb/No curves

  figure(1); 
  clf;
  semilogy(sim_in.EbNodB, fsk_theory.BERvec,'r;FSK theory;')
  hold on;
  semilogy(sim_in.EbNodB, fsk_sim.BERvec,'g;FSK sim;')
  semilogy(sim_in.EbNodB, fsk_fm_sim.BERvec,'b;FSK over FM sim;')
  hold off;
  grid("minor");
  axis([min(sim_in.EbNodB) max(sim_in.EbNodB) 1E-4 1])
  legend("boxoff");
  xlabel("Eb/No (dB)");
  ylabel("Bit Error Rate (BER)")

  % BER v SNR (3000 Hz noise BW and Eb=C/Rs=1/Rs)
  % Eb/No = (C/Rs)/(1/(N/B))
  % C/N   = (Eb/No)*(Rs/B)

  RsOnB_dB = 10*log10(sim_in.Rs/3000);
  figure(2); 
  clf;
  semilogy(sim_in.EbNodB+RsOnB_dB, fsk_theory.BERvec,'r;FSK theory;')
  hold on;
  semilogy(sim_in.EbNodB+RsOnB_dB, fsk_sim.BERvec,'g;FSK sim;')
  semilogy(sim_in.EbNodB+RsOnB_dB, fsk_fm_sim.BERvec,'b;FSK over FM sim;')
  hold off;
  grid("minor");
  axis([min(sim_in.EbNodB) max(sim_in.EbNodB) 1E-4 1])
  legend("boxoff");
  xlabel("S/N for RS=1200 bit/s and 3000 Hz noise bandwidth(dB)");
  ylabel("Bit Error Rate (BER)")
end

function run_fm_curves
  sim_in.nsam    = 96000;
  sim_in.verbose = 1;
  sim_in.pre_emp = 0;
  sim_in.de_emp  = 0;
  sim_in.CNdB = -4:2:20;

  sim_out = analog_fm_test(sim_in);

  figure(1)
  clf
  plot(sim_in.CNdB, sim_out.snrdB,"r;FM Simulated;");
  hold on;
  plot(sim_in.CNdB, sim_out.snr_theorydB,"g;FM Theory;");
  plot(sim_in.CNdB, sim_in.CNdB,"b; SSB Theory;");
  hold off;
  grid("minor");
  xlabel("FM demod input C/N (dB)");
  ylabel("FM demod output S/N (dB)");
  legend("boxoff");

  % C/No curves

  Bfm_dB = 10*log10(sim_out.Bfm);
  Bssb_dB = 10*log10(3000);

  figure(2)
  clf
  plot(sim_in.CNdB + Bfm_dB, sim_out.snrdB,"r;FM Simulated;");
  hold on;
  plot(sim_in.CNdB + Bfm_dB, sim_out.snr_theorydB,"g;FM Theory;");
  plot(sim_in.CNdB + Bssb_dB, sim_in.CNdB,"b; SSB Theory;");
  hold off;
  grid("minor");
  xlabel("FM demod input C/No (dB)");
  ylabel("FM demod output S/N (dB)");
  legend("boxoff");

endfunction

function run_fm_single
  sim_in.nsam    = 96000;
  sim_in.verbose = 2;
  sim_in.pre_emp = 0;
  sim_in.de_emp  = 0;

  sim_in.CNdB   = 0;
  sim_out = analog_fm_test(sim_in);
end

%run_fsk_curves
%run_fm_curves
run_fm_single

