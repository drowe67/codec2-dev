% tcohpsk.m
% David Rowe Oct 2014
%
% Octave script that tests the C port of the coherent PSK modem.  This
% script loads the output of unittest/tcohpsk.c and compares it to the
% output of the reference versions of the same functions written in
% Octave.

rand('state',1); 
randn('state',1);
graphics_toolkit ("gnuplot");

cohpsk;
autotest;

n = 160;
frames = 35;
framesize = 160;

load ../build_linux/unittest/tcohpsk_out.txt

sim_in = standard_init();
tx_bits = round(rand(1,framesize));
tx_bits_log = [];
for i=1:frames
  tx_bits_log = [tx_bits_log tx_bits];
end

stem_sig_and_error(1, 211, tx_bits_log_c(1:n), tx_bits_log(1:n) - tx_bits_log_c(1:n), 'tx bits', [1 n -1.5 1.5])
check(tx_bits_log, tx_bits_log_c, 'tx_bits');
