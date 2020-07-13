% fsk_lib_ldpc_tx.m
%
% LDPC coded 4FSK modem tx, generates a 8 kHz 16 bit short real valued sample file

function fsk_lib_ldpc_tx(filename, EbNodB=100, num_frames=10, Fs=8000, Rs=100)
  fsk_lib_ldpc;
  
  % set up LDPC code
  init_cml('~/cml/');
  load H_256_512_4.mat;
  [states code_param] = fsk_lib_ldpc_init (H, Rs, Fs);
  
  rand('seed',1);
  data_bits = round(rand(1,code_param.data_bits_per_frame)); tx_bits = [];
  for f=1:num_frames
    codeword_bits = LdpcEncode(data_bits, code_param.H_rows, code_param.P_matrix);
    tx_bits = [tx_bits codeword_bits];
  end

  tx = fsk_mod(states, tx_bits);

  % set up (optional) AWGN noise
  EcNodB = EbNodB + 10*log10(states.rate);
  EcNo = 10^(EcNodB/10);
  variance = states.Fs/(states.Rs*EcNo*states.bitspersymbol);

  % note real noise
  noise = sqrt(variance/2)*randn(length(tx),1);
  rx = tx + noise;

  frx=fopen(filename,"wb"); fwrite(frx, states.amp_scale*rx, "short"); fclose(frx);
  printf("Fs: %d Rs: %d EbNodB: %3.1f  EcNodB: %3.1f frames transmitted: %3d\n", Fs, Rs, EbNodB, EcNodB, num_frames);
end



