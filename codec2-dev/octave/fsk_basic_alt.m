% fsk_basic_alt.m
% David Rowe April 2018
%
% Modification of fsk_basic.m, to explore idea of alternating tone frequencies
% for multipath channels.

#{
  [ ] IL curves for half tone spacing
      + cf to channel with multipath
#}

% make sure we get same results from random number generators on each run

rand('seed',1);
randn('seed',1);

% constants -----------------------------------------------------------

Fs = 8000;     % sample rate
f1 = 1000;
f2 = 1800;
f3 = 1400;
f4 = 2200;
Rs = 800;      % symbol rate
Ts = Fs/Rs;    % length of each symbol in samples

% simulation parameters -----------------------------------------------

Nbits = 10000;
EbNodB = 9;    

% start simulation ----------------------------------------------------

tx_bits = round(rand(1,Nbits));

% continuous phase FSK modulator

w1 = 2*pi*f1/Fs;
w2 = 2*pi*f2/Fs;
w3 = 2*pi*f3/Fs;
w4 = 2*pi*f4/Fs;
tx_phase = 0;
tx = zeros(1,Ts*Nbits);

for i=1:Nbits
  for k=1:Ts
    if mod(i,2)
      if tx_bits(i)
        tx_phase += w2;
      else
        tx_phase += w1;
      end
    else
      if tx_bits(i)
        tx_phase += w4;
      else
        tx_phase += w3;
      end
    end    
    tx((i-1)*Ts+k) = exp(j*tx_phase);
  end
end

% AWGN channel noise

EbNo = 10^(EbNodB/10);
variance = Fs/(Rs*EbNo);
noise = sqrt(variance/2)*(randn(1,Nbits*Ts) + j*randn(1,Nbits*Ts));

d = Ts/2; l = length(tx);
rx = tx + noise;

figure(1); clf;
Tx = abs(fft(tx(1:Fs)));
plot(Tx);

% integrate and dump demodulator

rx_bits = zeros(1,Nbits);
for i=1:Nbits
  arx_symb = rx((i-1)*Ts + (1:Ts));
  if mod(i,2)
    filt1 = sum(exp(-j*w1*(1:Ts)) .* arx_symb);
    filt2 = sum(exp(-j*w2*(1:Ts)) .* arx_symb);
  else
    filt1 = sum(exp(-j*w3*(1:Ts)) .* arx_symb);
    filt2 = sum(exp(-j*w4*(1:Ts)) .* arx_symb);
  end
  rx_bits(i) = filt2 > filt1;
end

Nerrors = sum(xor(tx_bits, rx_bits));
ber = Nerrors/Nbits;
printf("EbNodB: %4.1f  Nerrors: %d BER: %1.3f\n", EbNodB, Nerrors, ber);

