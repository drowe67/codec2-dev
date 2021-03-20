% tofdm_acq.m
% Octave <-> C test for OFDM modem acquisition

ofdm_lib;
autotest;
randn('seed',1);
pkg load signal;

% generate a file of transmit samples
filename = "test_datac0.raw";
ofdm_tx(filename,"datac0",1,10,"awgn","bursts",1);

printf("\nRunning Octave version....\n");
config = ofdm_init_mode("datac0");
states = ofdm_init(config);  
states.verbose = 1; states.data_mode = "burst"; states.postambledectoren = 0;
ofdm_load_const;
frx=fopen(filename,"rb");
nin = states.nin;
rx = fread(frx, nin, "short")/(states.amp_scale/2);
f = 0;
while(length(rx) == nin)
  printf(" %2d ",f++);
  [timing_valid states] = ofdm_sync_search(states, rx);
  states.nin = nin;
  rx = fread(frx, nin, "short")/(states.amp_scale/2);
  printf("\n");
end   
fclose(frx);

printf("\nRunning C version....\n");
path_to_unittest = "../build_linux/unittest"
if getenv("PATH_TO_UNITEST")
  path_to_unittest_exe = getenv("PATH_TO_UNITTEST")
  printf("setting path from env var to %s\n", path_to_unittest);
end
system(sprintf("%s/tofdm_acq %s", path_to_unittest, filename));
load tofdm_acq_out.txt;

fg = 1;

% due to the order of processing we need to massage the Octave version a little before
% comparing
%tx_preamble = ofdm_clip(states, states.tx_preamble*states.amp_scale, states.ofdm_peak);
tx_preamble = ofdm_clip(states, states.tx_preamble, states.ofdm_peak);
stem_sig_and_error(fg, 211, real(tx_preamble_c), real(tx_preamble_c - tx_preamble), 'tx preamble re')
stem_sig_and_error(fg++, 212, imag(tx_preamble_c), imag(tx_preamble_c - tx_preamble), 'tx preamble im')
check(tx_preamble, tx_preamble_c, 'tx preamble', 0.1);


