% fsk_cont_phase.m
% David Rowe Feb 2016
%
% Looking at tx spectrum of FSK signal with and without continuous phase in the
% modulator.

1;

function tx  = fsk_mod(states, tx_bits)

    M  = states.M;
    Ts = states.Ts;
    Fs = states.Fs;
    ftx = states.ftx;
    phase_cont = states.phase_cont;

    num_bits = length(tx_bits);
    num_symbols = num_bits;
    tx = zeros(states.Ts*num_symbols,1);
    tx_phase = 0;
    s = 1;

    for i=1:num_bits

      % map bits to tone number

      tone = tx_bits(i) + 1;
 
      if phase_cont
        tx_phase_vec = tx_phase + (1:Ts)*2*pi*ftx(tone)/Fs;
        tx_phase = tx_phase_vec(Ts) - floor(tx_phase_vec(Ts)/(2*pi))*2*pi;
      else
        tx_phase_vec = (1:Ts)*2*pi*ftx(tone)/Fs;
      end

      tx((s-1)*Ts+1:s*Ts) = 2.0*cos(tx_phase_vec);
      s++;
      %printf("phase_cont: %d tx_phase_vec(Ts): %f\n", phase_cont, tx_phase_vec(Ts));

    end
endfunction


states.Fs  = 8000;
states.Rs  = 100;
states.Ts  = states.Fs/states.Rs;
states.M   = 2;
states.ftx = [1200 1305];  % need to choose these carefullto get disc phase

Nbits = 1024;

tx_bits = rand(1,Nbits) > 0.5;
%tx_bits = [0 1 0];

states.phase_cont = 1;
tx_cont  = fsk_mod(states, tx_bits);

states.phase_cont = 0;
tx_disc  = fsk_mod(states, tx_bits);

figure(1)
clf
plot(tx_cont(1:states.Ts*4))
hold on;
plot(tx_disc(1:states.Ts*4),'g')
hold off;
title('Time Domain')

figure(2)
clf
subplot(211)
Tx_cont = fft(tx_cont);
Tx_cont_dB = 20*log10(abs(Tx_cont));
plot(Tx_cont_dB)
axis([1 length(Tx_cont_dB)/2 0 100]);
grid
title('Cont Phase Spectra')

Tx_disc = fft(tx_disc);
Tx_disc_dB = 20*log10(abs(Tx_disc));
subplot(212)
plot(Tx_disc_dB,'g')
axis([1 length(Tx_cont_dB)/2 0 100]);
grid
title('Disc Cont Phase Spectra')
