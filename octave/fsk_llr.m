% fsk_llr.m
%
% 4FSK simulation to develop LLR estimation

#{
  TODO
  [X] Hard decision 4FSK modem simulation
  [X] single point AWGn simulation
  [ ] SD outputs
  [ ] histogram
  [ ] modify for frames of LDPC codec size
  [ ] LDPC codec with naive SD
#}

1;

function [rx_filt rx_bits] = run_single(tx_bits, M=4, EbNodB=100)
  bps = log2(M);  % bits per symbol
  Ts = 16;        % length of each symbol in samples

  Nbits = length(tx_bits);
  Nsymbols = Nbits/log2(M);

  mapper = bps:-1:1;
  % look up table demapper from symbols to bits (hard decision) 
  demapper=zeros(M,bps);
  for m=1:M
    for b=1:bps
      if  bitand(m-1,b) demapper(m,bps-b+1) = 1; end
    end
  end
  
  % continuous phase mFSK modulator

  w(1:M) = 2*pi*(1:M)/Ts;
  tx_phase = 0;
  tx = zeros(1,Ts*Nsymbols);

  for s=1:Nsymbols
    bits_for_this_symbol = tx_bits(bps*(s-1)+1:bps*s);
    symbol_index = bits_for_this_symbol * mapper' + 1;
    assert(demapper(symbol_index,:) == bits_for_this_symbol);
    for k=1:Ts
      tx_phase += w(symbol_index);
      tx((s-1)*Ts+k) = exp(j*tx_phase);
    end
  end

  % AWGN channel noise

  EsNodB = EbNodB + 10*log10(bps);
  EsNo = 10^(EsNodB/10);
  variance = Ts/EsNo;
  noise = sqrt(variance/2)*(randn(1,Nsymbols*Ts) + j*randn(1,Nsymbols*Ts));
  rx = tx + noise;

  % integrate and dump demodulator

  rx_bits = zeros(1,Nbits);
  rx_filt = zeros(Nsymbols,M);
  for s=1:Nsymbols
    arx_symb = rx((s-1)*Ts + (1:Ts));
    for m=1:M
      rx_filt(s,m) = abs(sum(exp(-j*w(m)*(1:Ts)) .* arx_symb));
    end
    [tmp symbol_index] = max(rx_filt(s,:));
    rx_bits(bps*(s-1)+1:bps*s) = demapper(symbol_index,:);
  end

  Nerrors = sum(xor(tx_bits, rx_bits));
  ber = Nerrors/Nbits;
  printf("EbNodB: %4.1f  Nerrors: %d BER: %1.3f\n", EbNodB, Nerrors, ber);
endfunction

function plot_hist(Nbits,EbNodB)
  Nbits = 10000; M=4;
  tx_bits = ones(1,Nbits);
  rx_filt = run_single(tx_bits,M,EbNodB=6);
  figure(1); clf;
  for m=1:4
    subplot(2,2,m);
    hist(rx_filt(:,m),25);
  end
endfunction

rand('seed',1);
randn('seed',1);
format short
more off

% Eb/No = 6dB test pount should be about BER = 0.0157
Nbits = 10000; tx_bits = round(rand(1,Nbits)); run_single(tx_bits,4,6);

plot_hist;
