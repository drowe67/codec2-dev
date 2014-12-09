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

rand('state',1); 
randn('state',1);

function sim_out = ber_test(sim_in)
  Fs        = 96000;
  fmark     = 1200;
  fspace    = 2200;
  Rs        = sim_in.Rs;
  Ts        = Fs/Rs;

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
  Fs        = 96000;
  fm_max    = 3000;
  fm        = 1000; wm = 2*pi*fm/Fs;

  nsam      = sim_in.nsam;
  CNdB      = sim_in.CNdB;

  % FM modulator constants

  fc = 24E3; wc = 2*pi*fc/Fs;
  fd = 3E3; wd = 2*pi*fd/Fs;
  Bfm = 2*(fd+fm_max);             % Carson's rule for FM demod input noise BW
  printf("Bfm: %f\n", Bfm);
  bin  = fir1(100,(Bfm/2)/(Fs/2));
  bout = fir1(100,fm_max/(Fs/2));
  % simulate over a range of C/N values

  for ne = 1:length(CNdB)

    aCNdB = CNdB(ne);
    CN = 10^(aCNdB/10);
    variance = Fs/(CN*Bfm);

    % FM Modulator -------------------------------

    t = 0:(nsam-1);
    tx = exp(j*(wc*t + wd*cos(wm*t)));

    % Channel ---------------------------------

    noise = sqrt(variance/2)*(randn(1,length(tx)) + j*randn(1,length(tx)));
    rx    = (1 + 0.0*cos(2*pi*1000*t/Fs)) .* tx + noise;
    %printf("p rx: %f\n", var(rx))
    figure(3)
    subplot(211)
    plot(20*log10(abs(fft(rx))))
    title('FM Demod input Spectrum');
    axis([1 length(tx) 0 100]);
    
    % FM Demodulator

    l = length(rx);
    rx_bb = rx .* exp(-j*wc*t);               % down to complex baseband
    %printf("p rx_bb: %f\n", var(rx_bb))
    %p1 = (rx_bb * rx_bb')/l;

    figure(1)
    subplot(211)
    plot(rx_bb(100:1000),'+');
    title('FM demod input')
    axis([-2 2 -2 2]);

    rx_bb = filter(bin,1,rx_bb);
    %rx_bb = rx_bb .* (1.0+0.2*randn(1,l));
    p2 = (rx_bb * rx_bb')/l;
    %printf("p rx_bb filter: %f C/N: %f\n C/N dB: %f", p2, 1/p2, 10*log10(1/p2))

    subplot(212)
    plot(rx_bb(100:1000),'+');
    axis([-2 2 -2 2]);
 
    figure(3)
    subplot(212)
    Rx_bb = 20*log10(abs(fft(rx_bb)));
    plot(Rx_bb)
    title('FM Demod input Spectrum filter');

    %rx_bb_diff = [ 1 rx_bb(2:l) .* conj(rx_bb(1:l-1))];
    rx = atan2(imag(rx_bb),real(rx_bb));
    rx = filter(bout,1,rx);
    w = 2*pi*fm/Fs; beta = 0.99;
    rx = filter([1 -2 1],[1 -2*beta beta*beta], rx);

    % notch out test tone

    w = 2*pi*fm/Fs; beta = 0.99;
    rx_notch = filter([1 -2*cos(w) 1],[1 -2*beta*cos(w) beta*beta], rx);

    figure(2)
    subplot(211)
    plot(rx_notch)
    %axis([1 800 -0.3 0.3]) 
    subplot(212)
    Rx = 20*log10(abs(fft(rx_notch(1000:l))));
    plot(Rx(1:1100))
    p3 = (rx(1000:l) * rx(1000:l)')/(l-1000);
    p4 = (rx_notch(1000:l) * rx_notch(1000:l)')/(l-1000);

    %printf("%f %f\n", aCNdB, CN);
    %sig = sqrt(2)*cos(2*pi*wm*t) + sqrt(1/CN)*randn(1,l);
    %p3  = (sig * sig')/l;
    snr = (p3-p4)/p4;
    printf("C/N: %f S+N: %f  N: %f SNR: %f\n", aCNdB, p3, p4, 10*log10(snr));
    sim_out.p3(ne) = p3;
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

sim_in.nsam = 96000;
sim_in.CNdB = [20 30 40 50 100];
sim_out = analog_fm_test(sim_in);

if 0
%snrdB_in = 10:2:30
snrdB_in = 30;
%sim_in.CNdB = [200 snrdB_in];

s  = sim_out.p3(1)
sn = sim_out.p3(2:length(sim_out.p3))
n  = sn - s;
snrdB_out = 10*log10(s ./ n)

figure(5)
plot(snrdB_in, snrdB_out)
end
