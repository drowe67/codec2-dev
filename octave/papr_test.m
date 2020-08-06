% papr_test.m
%
#{

  TODO:
    [ ] measure BER
    [ ] curves with different clipping
    [ ] companding
    [ ] heat map type scatter diagram
#}

M = 160;      % length of one symbol in samples
Nc = 8;      % number of carriers
Nsym = 1000;   % number of symbols to simulate
bps  = 2;     % two bits per symbol for QPSK

function symbol = qpsk_mod(two_bits)
    two_bits_decimal = sum(two_bits .* [2 1]); 
    switch(two_bits_decimal)
        case (0) symbol =  1;
        case (1) symbol =  j;
        case (2) symbol = -j;
        case (3) symbol = -1;
    endswitch
endfunction

function two_bits = qpsk_demod(symbol)
    bit0 = real(symbol*exp(j*pi/4)) < 0;
    bit1 = imag(symbol*exp(j*pi/4)) < 0;
    two_bits = [bit1 bit0];
endfunction

% generate a 2D array of QPSK symbols

Nphases = 2^bps;
tx_phases = pi/2*floor((rand(Nsym,Nc)*Nphases));
tx_sym = exp(j*tx_phases);

w = 2*pi/M;
tx = [];

% generate OFDM signal

for s=1:Nsym
  atx = zeros(1,M);
  for c=1:Nc
    atx += exp(j*(0:M-1)*c*w)*tx_sym(s,c);
  end
  tx = [tx atx];
end
Nsam = length(tx);

% AWGN channel

EbNodB = 4;
EsNodB = EbNodB + 10*log10(bps);
variance = M/(10^(EsNodB/10));
noise = sqrt(variance/2)*randn(1,Nsam) + j*sqrt(variance/2)*randn(1,Nsam);
rx = tx + noise;

% demodulate

rx_sym = zeros(Nsym,Nc);
for s=1:Nsym
  st = (s-1)*M+1; en = s*M;
  for c=1:Nc
    rx_sym(s,c) = sum(exp(-j*(0:M-1)*c*w) .* rx(st:en))/M;
  end
end

% count bit errors

Tbits = Terrs = 0;
for s=1:Nsym
  for c=1:Nc
    tx_bits = qpsk_demod(tx_sym(s,c));
    rx_bits = qpsk_demod(rx_sym(s,c));
    Tbits += bps;
    Terrs += sum(xor(tx_bits,rx_bits));
  end
end

figure(1); clf;
plot(abs(tx(1:10*M)))
figure(2); clf; hist(abs(tx),25)
figure(3); clf; plot(rx_sym,'+'); axis([-2 2 -2 2]);
papr = 20*log10(max(abs(tx))/mean(abs(tx)));
printf("PAPR: %5.2f dB Tbits: %d Terrs: %d BER: %5.3f\n", papr, Tbits, Terrs, Terrs/Tbits)
  
