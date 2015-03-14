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
sim_in.framesize = 160;
sim_in.ldpc_code = 0;
sim_in.ldpc_code_rate = 1;
sim_in.Nc               = 4;
sim_in.Rs               = 50;
sim_in.Ns               = 4;
sim_in.Np               = 2;
sim_in.Nchip            = 1;
sim_in.modulation       = 'qpsk';
sim_in.do_write_pilot_file = 0;
sim_in = symbol_rate_init(sim_in);

rand('state',1); 
tx_bits = round(rand(1,framesize));
tx_bits_log = [];
tx_symb_log = [];
for i=1:frames
  tx_bits_log = [tx_bits_log tx_bits];
  [tx_symb tx_bits prev_tx_sym] = bits_to_qpsk_symbols(sim_in, tx_bits, [], []);
  tx_symb_log = [tx_symb_log; tx_symb];
end

stem_sig_and_error(1, 111, tx_bits_log_c(1:n), tx_bits_log(1:n) - tx_bits_log_c(1:n), 'tx bits', [1 n -1.5 1.5])
stem_sig_and_error(2, 211, real(tx_symb_log_c(1:n)), real(tx_symb_log(1:n) - tx_symb_log_c(1:n)), 'tx symb re', [1 n -1.5 1.5])
stem_sig_and_error(2, 212, imag(tx_symb_log_c(1:n)), imag(tx_symb_log(1:n) - tx_symb_log_c(1:n)), 'tx symb im', [1 n -1.5 1.5])
check(tx_bits_log, tx_bits_log_c, 'tx_bits');
check(tx_symb_log, tx_symb_log_c, 'tx_symb');
