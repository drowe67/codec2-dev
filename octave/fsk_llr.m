% fsk_llr.m
%
% 4FSK simulation to develop LLR estimation

#{
  TODO
  [ ] Hard decision 4FSK modme simulation
  [ ] single point AWGn simulation
  [ ] SD outputs
  [ ] histogram
  [ ] modify for frames of LDPC codec size
  [ ] LDPC codec with naive SD
#}

1;

function run_single(EbNodB)
  rand('seed',1);
  randn('seed',1);
  format short
  more off

  M = 4;          % M-FSK
  bps = log2(M);  % bits per symbol
  Ts = 16;        % length of each symbol in samples
  Nbits = 10;

  tx_bits = round(rand(1,Nbits));
  Nsymbols = Nbits/log2(M);

  symbol_lut=zeros(M,bps);
  for m=1:M
    for b=1:bps
      if  bitand(m-1,b) symbol_lut(m,bps-b+1) = 1; end
    end
  end

  % continuous phase mFSK modulator

  w(1:M) = 2*pi*(1:M)/Ts;
  tx_phase = 0;
  tx = zeros(1,Ts*Nsymbols);

  for s=1:Nsymbols
    bits_for_this_symbol = tx_bits(bps*(s-1)+1:bps*s);
    symbol_index = bits_for_this_symbol * [2 1]' + 1;
    assert(symbol_lut(symbol_index,:) == bits_for_this_symbol);
    for k=1:Ts
      tx_phase += w(symbol_index);
      tx((s-1)*Ts+k) = exp(j*tx_phase);
    end
  end

  % AWGN channel noise

  EbNo = 10^(EbNodB/10);
  variance = Ts/(EbNo*bps);
  noise = sqrt(variance/2)*(randn(1,Nsymbols*Ts) + j*randn(1,Nsymbols*Ts));
  rx = tx + noise;

  % integrate and dump demodulator

  rx_bits = zeros(1,Nbits);
  filt = zeros(1,M);
  for s=1:Nsymbols
    arx_symb = rx((s-1)*Ts + (1:Ts));
    for m=1:M
      filt(m) = sum(exp(-j*w(m)*(1:Ts)) .* arx_symb);
    end
    [tmp symbol_index] = max(filt);
    rx_bits(bps*(s-1)+1:bps*s) = symbol_lut(symbol_index,:);
  end

  Nerrors = sum(xor(tx_bits, rx_bits));
  ber = Nerrors/Nbits;
  printf("EbNodB: %4.1f  Nerrors: %d BER: %1.3f\n", EbNodB, Nerrors, ber);
endfunction

run_single(100)
