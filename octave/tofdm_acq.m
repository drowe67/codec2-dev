% tofdm_acq.m
% Octave <-> C test for OFDM modem acquisition

ofdm_lib;
autotest;
randn('seed',1);
pkg load signal;
more off;

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
timing_mx_log = []; ct_est_log = []; foff_est_log = []; timing_valid_log = [];
while(length(rx) == nin)
  printf(" %2d ",f++);
  [timing_valid states] = ofdm_sync_search(states, rx);
  states.nin = nin;
  timing_mx_log = [timing_mx_log states.timing_mx];
  ct_est_log = [ct_est_log states.ct_est];
  foff_est_log = [foff_est_log states.foff_est_hz];
  timing_valid_log = [timing_valid_log states.timing_valid];
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

fg = 1; p = 0;

tx_preamble = states.tx_preamble;
stem_sig_and_error(fg, 211, real(tx_preamble_c), real(tx_preamble_c - tx_preamble), 'tx preamble re')
stem_sig_and_error(fg++, 212, imag(tx_preamble_c), imag(tx_preamble_c - tx_preamble), 'tx preamble im')
p += check(tx_preamble, tx_preamble_c, 'tx preamble', 0.1);

stem_sig_and_error(fg, 211, real(timing_mx_log_c), real(timing_mx_log_c - timing_mx_log), 'timing mx')
p += check(timing_mx_log, timing_mx_log_c, 'timing_mx');
stem_sig_and_error(fg++, 212, real(ct_est_log_c), real(ct_est_log_c - ct_est_log), 'ct est')
p += check(ct_est_log, ct_est_log_c, 'ct_est_mx');

stem_sig_and_error(fg, 211, real(foff_est_log_c), real(foff_est_log_c - foff_est_log), 'foff est')
p += check(foff_est_log, foff_est_log_c, 'foff_est');
stem_sig_and_error(fg++, 212, real(timing_valid_log_c), real(timing_valid_log_c - timing_valid_log), 'timing valid')
p += check(timing_valid_log, timing_valid_log_c, 'timing_valid');

if p == 5 printf("PASS\n"); else printf("FAIL\n"); end

  



