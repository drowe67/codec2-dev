% tofdm.m
% David Rowe and Steve Sampson June 2017
%
% Octave script for comparing Octave and C versions of OFDZM modem

Frames = 10;

more off;
ofdm_lib;
autotest;

% Run a few frames of Octave version

Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 16; Ns = 8;
states = ofdm_initbps, Rs, Tcp, Ns, Nc);
ofdm_load_const;

rand('seed',1);
tx_bits = round(rand(1,Nbitsperframe));

tx_bits_log = []; tx_log = [];
for f=1:Frames
  tx_bits_log = [tx_bits tx_bits];
  tx_log = [tx_log ofdm_mod(states, tx_bits)];
end
 
% Load C version and plot Octave and C states and differences -----------------------------

load ../build_linux/octave/tofdm_out.txt;

stem_sig_and_error(1, 111, tx_bits_log_c, tx_bits_log - tx_bits_log_c, 'tx bits', [1 length(tx_bits_log) -1.5 1.5])
stem_sig_and_error(2, 211, real(tx_log_c), real(tx_log - tx_log_c), 'tx re', [1 length(tx_log_c) -1.5 1.5])
stem_sig_and_error(2, 212, imag(tx_log_c), real(tx_log - tx_log_c), 'tx im', [1 length(tx_log_c) -1.5 1.5])

% Run through checklist -----------------------------

check(W, W_c, 'W');
check(tx_bits_log, tx_bits_log_c, 'tx_bits');
check(tx_symbols_log, tx_log_c, 'tx');
