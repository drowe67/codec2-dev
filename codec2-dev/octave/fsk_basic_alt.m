% fsk_basic_alt.m
% David Rowe April 2018
%
% Modification of fsk_basic.m, to explore idea of alternating tone frequencies
% for multipath channels.

#{

  So alternate tone spacings at Rs/2 worked quite well, no implementation loss
  but odd spectrum.

  But when multipath added (e.g 2ms delay on 2.5ms Rs), performance the same
  as for non-stepped tones.  Multipath echo must be interfering with stepped
  tone.

  When 2nd step of tones are Rs spaced - then multipath perf much better,
  as expected.

  There is also the impact of freq selective fading as well as ISI, FSK
  tone level will be pushed down and we may choose wrong tone.

  Also be interesting to test FSK with FEC on HF, combined response.
#}

% make sure we get same results from random number generators on each run

rand('seed',1);
randn('seed',1);

% constants -----------------------------------------------------------

Fs = 8000;     % sample rate
f1 = 1000;
f2 = 1400;
f3 = 1200;
f4 = 1600;
Rs = 400;      % symbol rate
Ts = Fs/Rs;    % length of each symbol in samples

% simulation parameters -----------------------------------------------

Nbits = 10000;
EbNodB = 9;    
multipath = 1;
stepped = 1;

% start simulation ----------------------------------------------------

tx_bits = round(rand(1,Nbits));

% continuous phase FSK modulator

w1 = 2*pi*f1/Fs;
w2 = 2*pi*f2/Fs;
w3 = 2*pi*f3/Fs;
w4 = 2*pi*f4/Fs;

% choose tones here, e.g. offset on every 2nd symbol

if stepped
  w = [w1 w2; w3 w4];
else
  w = [w1 w2; w1 w2];
end

tx_phase = 0;
tx = zeros(1,Ts*Nbits);

for i=1:Nbits
  for k=1:Ts
    tx_phase += w(1+mod(i-1,2),1+tx_bits(i));
    tx((i-1)*Ts+k) = exp(j*tx_phase);
  end
end

% AWGN channel noise

EbNo = 10^(EbNodB/10);
variance = Fs/(Rs*EbNo);
noise = sqrt(variance/2)*(randn(1,Nbits*Ts) + j*randn(1,Nbits*Ts));

rx = tx + noise;
if multipath
  d = 0.002*Fs;
  tx_delayed = [zeros(1,d) tx(1:end-d)];
  rx +=  0.5*tx_delayed;
end

figure(1); clf;
Tx = abs(fft(tx(1:Fs)));
plot(Tx);

% integrate and dump demodulator

rx_bits = zeros(1,Nbits);
for i=1:Nbits
  arx_symb = rx((i-1)*Ts + (1:Ts));
  wr = 1+mod(i-1,2);
  filt1 = sum(exp(-j*w(wr,1)*(1:Ts)) .* arx_symb);
  filt2 = sum(exp(-j*w(wr,2)*(1:Ts)) .* arx_symb);
  rx_bits(i) = filt2 > filt1;
end

Nerrors = sum(xor(tx_bits, rx_bits));
ber = Nerrors/Nbits;
printf("EbNodB: %4.1f  Nerrors: %d BER: %1.3f\n", EbNodB, Nerrors, ber);

