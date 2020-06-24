% fsk_llr.m
%
% 4FSK simulation to develop LLR estimation algorithms for 4FSK/LDPC modems

#{
  TODO
  [X] Hard decision 4FSK modem simulation
  [X] single point AWGn simulation
  [X] filter outputs
  [X] histogram
  [X] modify for frames of LDPC codec size
  [X] LDPC codec with naive SD
#}

1;
ldpc;

% single Eb/No point simulation

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
  printf("EbNodB: %4.1f  Uncoded Nbits: %5d Nerrors: %4d BER: %1.3f\n", EbNodB, Nbits, Nerrors, ber);
endfunction


% Plot histograms of Rx filter outputs

function plot_hist(Nbits=10000,EbNodB=2)
  M=4;
  tx_bits = ones(1,Nbits);
  rx_filt = run_single(tx_bits,M,EbNodB);
  figure(1); clf; X = 0:31;
  for m=1:4
    subplot(2,2,m);
    H=hist(rx_filt(:,m),X);
    plot(X,H);
  end
endfunction


% 2FSK SD -> LLR mapping that we used for Wenet SSTV system

function llr = sd_to_llr(sd)
    sd = sd / mean(abs(sd));
    x = sd - sign(sd);
    sumsq = sum(x.^2);
    summ = sum(x);
    mn = summ/length(sd);
    estvar = sumsq/length(sd) - mn*mn; 
    estEsN0 = 1/(2* estvar + 1E-3); 
    llr = 4 * estEsN0 * sd;
endfunction


% single point LDPC encoded frame simulation, usin 2FSK as a tractable starting point
% Note: ~/cml/matCreateConstellation.m has some support for FSK - can it do 4FSK?

function run_single_ldpc(Nbits=10000,EbNodB=4)
  M=2;
  bps = 1; modulation = 'FSK'; mod_order=2; mapping = 'gray'; decoder_type = 0; max_iterations = 100;
  load H_256_768_22.txt
  code_param = ldpc_init_user(H_256_768_22, modulation, mod_order, mapping);
  Nframes = floor(Nbits/code_param.data_bits_per_frame);
  Nbits = Nframes*code_param.data_bits_per_frame;

  % Encoder
  data_bits = round(rand(1,code_param.data_bits_per_frame));
  tx_bits = [];
  for f=1:Nframes;
    codeword_bits = LdpcEncode(data_bits, code_param.H_rows, code_param.P_matrix);
    tx_bits = [tx_bits codeword_bits];
  end  

  % modem/channel simulation
  rx_filt = run_single(tx_bits,M,EbNodB);

  % Decoder
  Nerrors = 0;
  for f=1:Nframes
    st = (f-1)*code_param.coded_bits_per_frame + 1;
    en = st + code_param.coded_bits_per_frame - 1;
    sd = rx_filt(st:en,1) - rx_filt(st:en,2);
    llr = sd_to_llr(sd)';
  
    [x_hat, PCcnt] = MpDecode(llr, code_param.H_rows, code_param.H_cols, ...
                            max_iterations, decoder_type, 1, 1);         
    Niters = sum(PCcnt!=0);
    detected_data = x_hat(Niters,:);
    Nerrors += sum(xor(data_bits, detected_data(1:code_param.data_bits_per_frame)));
  end  
  ber = Nerrors/Nbits;
  printf("EbNodB: %4.1f  Coded   Nbits: %5d Nerrors: %4d BER: %1.3f\n", EbNodB, Nbits, Nerrors, ber);
endfunction

% Choose what you would like to run here --------------------------------------------------

rand('seed',1);
randn('seed',1);
format short
more off
init_cml('~/cml/');

% 1) Eb/No = 6dB test pount should be about BER = 0.0157
Nbits = 10000; tx_bits = round(rand(1,Nbits)); run_single(tx_bits,4,6);

% 2) Histograms
plot_hist;

% 3) 2FSK/LDPC simulation as a starting point
run_single_ldpc;

