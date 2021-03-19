% tofdm_acq.m
% Octave <-> C test for OFDM moedm acquisition

ofdm_lib;
autotest;
randn('seed',1);
pkg load signal;

config = ofdm_init_mode("datac0");
states = ofdm_init(config);  
print_config(states);
ofdm_load_const;

printf("\nRunning C version....\n");
path_to_unittest = "../build_linux/unittest"
if getenv("PATH_TO_UNITEST")
  path_to_unittest_exe = getenv("PATH_TO_UNITTEST")
  printf("setting path from env var to %s\n", path_to_unittest);
end
system(sprintf("%s/tofdm_acq", path_to_unittest));
load tofdm_acq_out.txt;

fg = 1;

% due to the order of processing we need to massage the Octave version a little before
% comparing
tx_preamble = ofdm_clip(states, states.tx_preamble*states.amp_scale, states.ofdm_peak);
stem_sig_and_error(fg, 211, real(tx_preamble_c), real(tx_preamble_c - tx_preamble), 'tx preamble re')
stem_sig_and_error(fg++, 212, imag(tx_preamble_c), imag(tx_preamble_c - tx_preamble), 'tx preamble im')
check(tx_preamble, tx_preamble_c, 'tx preamble', 0.1);


