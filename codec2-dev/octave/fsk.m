% fsk.m
% David Rowe Nov 2014

% Simulation to test FSK demod
%
% TODO
%   [X] Code up mod/non-coh demod/AWGN channel simulation
%   [X] Eb/No verses BER curves
%   [ ] test with pre/de-emphahsis impairments
%       + this will introduce delay, use fir filter, group delay
%   [ ] test with hard limiting      
%       + it's FSK, so AM shouldn't matter?  
%       + also make 8-bit fixed point impl easy
%   [ ] channel simulation of HT/FM radio
%       + filtering, varying modulation index
%   [ ] GMSK
%   [ ] refactor to plot analog FM demod curves
%   [ ] SSB curves
%   [ ] different modn index beta curves, 
%   [ ]  illustration of harmonic dist, 
%   [ ] integration with AFSK, AFSK-FM, v AFSK-SSB (ie FSK)


rand('state',1); 
randn('state',1);
graphics_toolkit ("gnuplot");

function sim_out = ber_test(sim_in)
  Fs        = 96000;
  fmark     = 1200;
  fspace    = 2200;
  Rs        = sim_in.Rs;
  Ts        = Fs/Rs;
  emphasis  = 50E-6;

  framesize = sim_in.framesize;
  EbNodB    = sim_in.EbNodB;

  % FM modulator constants

  fc = 24E3; wc = 2*pi*fc/Fs;
  f_max_deviation = 3E3; w_max_deviation = 2*pi*f_max_deviation/Fs

  % simulate over a rnage of Eb/No values

  for ne = 1:length(EbNodB)
    Nerrs = Terrs = Tbits = 0;

    aEbNodB = EbNodB(ne);
    EbNo = 10^(aEbNodB/10);
    variance = Fs/(Rs*EbNo)

    % Modulator -------------------------------

    tx_bits = round(rand(1, framesize));
    %tx_bits = zeros(1, framesize);
    tx = zeros(1,framesize*Ts);
    tx_phase = 0;

    for i=1:framesize
      for k=1:Ts
        if tx_bits(i) == 1
          tx_phase += 2*pi*fmark/Fs;
        else
          tx_phase += 2*pi*fspace/Fs;
        end
        tx_phase = tx_phase - floor(tx_phase/(2*pi))*2*pi;
        tx((i-1)*Ts+k) = sqrt(2)*cos(tx_phase);
      end
    end

    figure(4)
    subplot(211);
    plot(tx(1:1000));

    % Optional AFSK-FM modulator

    if sim_in.fm
      mod = tx/sqrt(2);
      wt = wc*(0:(framesize*Ts-1)) + w_max_deviation*mod;
      % wt = wc*(0:(framesize*Ts-1));
      tx = exp(j*wt);
    end
  
    % Channel ---------------------------------

    noise = sqrt(variance/2)*(randn(1,length(tx)) + j*randn(1,length(tx)));
    rx    = tx + noise;
    printf("EbNo: %f Eb: %f var No: %f EbNo (meas): %f\n", 
    EbNo, var(tx)*Ts/Fs, var(noise)/(Fs/2), (var(tx)*Ts/Fs)/(var(noise)/(Fs/2)));

    %rx = exp(j*angle(rx));
    figure(4)
    subplot(211)
    plot(real(rx(1:1000)),'+');
    title('FM demod input (real)')

    % Optional AFSK-FM demodulator

    if sim_in.fm
      figure(4)
      subplot(211)
      plot(real(rx(1:1000)),'+');
      title('FM demod input (real)')

      l = length(rx);
      rx_bb = rx .* exp(-j*wc*(0:(l-1)));               % down to complex baseband
      b=fir1(15,f_max_deviation/Fs);
      rx_bb=filter(b,1,rx_bb);
      rx_bb_diff = rx_bb(2:l) .* conj(rx_bb(1:l-1));    % difference in phase, which is freq 
                                                        % also keeps us away from atan function 
                                                        % discontinuities

      rx = atan2(imag(rx_bb),real(rx_bb));
      subplot(212);
      plot(rx(1:1000));
      title('FM demod output');

      figure(5)
      plot(rx_bb)
      axis([-2 2 -2 2])   
    end

    % Demodulator -----------------------------

    % non-coherent FSK demod

    mark_dc  = rx .* exp(-j*(0:(framesize*Ts)-1)*2*pi*fmark/Fs);
    space_dc = rx .* exp(-j*(0:(framesize*Ts)-1)*2*pi*fspace/Fs);

    rx_bits = zeros(1, framesize);
    for i=1:framesize
      st = (i-1)*Ts+1;
      en = st+Ts-1;
      mark_int(i)  = sum(mark_dc(st:en));
      space_int(i) = sum(space_dc(st:en));
      rx_bits(i) = abs(mark_int(i)) > abs(space_int(i));
    end
  
    figure(2);
    subplot(211)
    stem(abs(mark_int(1:20)));
    subplot(212)
    stem(abs(space_int(1:20)));
    
    error_positions = xor(rx_bits, tx_bits);
    Nerrs = sum(error_positions);
    Terrs += Nerrs;
    Tbits += framesize-1;

    TERvec(ne) = Terrs;
    BERvec(ne) = Terrs/Tbits;

    printf("EbNo (db): %3.2f Terrs: %d BER: %3.2f \n", aEbNodB, Terrs, Terrs/Tbits);
  end

  sim_out.TERvec = TERvec;
  sim_out.BERvec = BERvec;
endfunction


function sim_out = analog_fm_test(sim_in)
  Fs        = 96000; FsOn2 = Fs/2;
  fm        = 1000; wm = 2*pi*fm/Fs;
  nsam      = Fs;
  CNdB      = sim_in.CNdB;
  verbose   = sim_in.verbose;
  pre_emp   = sim_in.pre_emp;
  de_emp    = sim_in.de_emp;

  % FM modulator constants

  fm_max = 3000;                   % max modulation freq
  fc = 24E3; wc = 2*pi*fc/Fs;      % carrier frequency
  fd = 5E3; wd = 2*pi*fd/Fs;       % (max) deviation
  m = fd/fm_max;                   % modulation index
  Bfm = 2*(fd+fm_max);             % Carson's rule for FM signal bandwidth
  tc  = 50E-6;
  prede = [1 -(1 - 1/(tc*Fs))];    % pre/de emp filter coeffs

  %printf("Bfm: %f m: %f wc: %f wd: %f wm: %f\n", Bfm, fd/fm_max, wc, wd, wm);
  
  % input filter gets rid of excess noise before demodulator, as too much
  % noise causes atan2() to jump around, e.g. -pi to pi.  However this
  % filter can cause harmonic distortion at very high SNRs, as it knocks out
  % some of the FM signal spectra.  This filter isn't really required for high
  % SNRs > 20dB.

  fc = (Bfm/2)/(FsOn2);
  bin  = firls(200,[0 fc*(1-0.05) fc*(1+0.05) 1],[1 1 0.01 0.01]);

  % demoduator output filter to limit us to 3 kHz

  fc = (fm_max)/(FsOn2);
  bout = firls(200,[0 0.95*fc 1.05*fc 1], [1 1 0.01 0.01]);

  % start simulation

  for ne = 1:length(CNdB)

    % work out the variance we need to obtain our C/N in the bandwidth
    % of the FM demod.  The gaussian generator randn() generates noise
    % with a bandwidth of Fs

    aCNdB = CNdB(ne);
    CN = 10^(aCNdB/10);
    variance = Fs/(CN*Bfm);
     
    % FM Modulator -------------------------------

    t = 0:(nsam-1);
    tx_phase = 0;
    mod = sin(wm*t);
    if pre_emp
      mod = filter(prede,1,mod);
      mod = mod/max(mod);           % AGC to set deviation
    end
    for i=0:nsam-1
        w = wc + wd*mod(i+1);
        tx_phase = tx_phase + w;
        tx_phase = tx_phase - floor(tx_phase/(2*pi))*2*pi;
        tx(i+1) = exp(j*tx_phase);
    end       

    % Channel ---------------------------------

    noise = sqrt(variance/2)*(randn(1,nsam) + j*randn(1,nsam));
    rx = tx + noise;

    % FM Demodulator

    rx_bb = rx .* exp(-j*wc*t);      % down to complex baseband
    rx_bb = filter(bin,1,rx_bb);
    p2 = (rx_bb * rx_bb')/nsam;
 
    rx_bb_diff = [ 1 rx_bb(2:nsam) .* conj(rx_bb(1:nsam-1))];
    rx_out = atan2(imag(rx_bb_diff),real(rx_bb_diff));
    rx_out = filter(bout,1,rx_out);
    if de_emp
      rx_out = filter(1,prede,rx_out);
    end

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

    %hist(rx_notch(1000:l),100)
  end

endfunction

more off;

function run_fsk_sim
  sim_in.Rs        = 1200;
  sim_in.framesize = 1200;
  sim_in.EbNodB    = 30;
  sim_in.fm        = 1;

  EbNo  = 10 .^ (sim_in.EbNodB/10);
  non_coh_theory.BERvec = 0.5*exp(-EbNo/2);
  non_coh_sim = ber_test(sim_in);

  figure(1); 
  clf;
  semilogy(sim_in.EbNodB, non_coh_theory.BERvec,'r;FSK non coherent AWGN theory;')
  hold on;
  semilogy(sim_in.EbNodB, non_coh_sim.BERvec,'g;FSK non coherent AWGN sim;')
  hold off;
  grid("minor");
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
  plot(sim_in.CNdB, sim_in.CNdB,"g; SSB Theory;");
  hold off;
  grid;
  xlabel("FM demod input CNR (dB)");
  ylabel("FM demod output SNR (dB)");
endfunction

function run_fm_single
  sim_in.nsam    = 96000;
  sim_in.verbose = 2;
  sim_in.pre_emp = 0;
  sim_in.de_emp  = 0;

  sim_in.CNdB   = 0;
  sim_out = analog_fm_test(sim_in);
end

%run_fm_curves
run_fm_single

