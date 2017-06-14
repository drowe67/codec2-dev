% tofdm.m
% David Rowe and Steve Sampson June 2017
%
% Octave script for comparing Octave and C versions of OFDZM modem


Nframes = 2;

more off;
ofdm_lib;
autotest;

% ---------------------------------------------------------------------
% Run Octave version 
% ---------------------------------------------------------------------

Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 16; Ns = 8;
states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
ofdm_load_const;

rand('seed',1);
tx_bits = round(rand(1,Nbitsperframe));

% Run tx loop

tx_bits_log = []; tx_log = [];
for f=1:Nframes
  tx_bits_log = [tx_bits_log tx_bits];
  tx_log = [tx_log ofdm_mod(states, tx_bits)];
end

% Channel simulation

rx = tx_log;

% Init rx with ideal timing so we can test with timing estimation disabled

Nsam = length(rx);
prx = 1;
nin = Nsamperframe+2*(M+Ncp);
states.rxbuf(Nrxbuf-nin+1:Nrxbuf) = rx(prx:nin);
prx += nin;

rxbuf_log = [];

for f=1:Nframes

  % insert samples at end of buffer, set to zero if no samples
  % available to disable phase estimation on future pilots on last
  % frame of simulation

  lnew = min(Nsam-prx,states.nin);
  rxbuf_in = zeros(1,states.nin);

  if lnew
    rxbuf_in(1:lnew) = rx(prx:prx+lnew-1);
  end
  prx += states.nin;
  [rx_bits_raw states aphase_est_pilot_log arx_np arx_amp] = ofdm_demod(states, rxbuf_in);

  % log some states for comparison to C

  rxbuf_log = [rxbuf_log states.rxbuf];
end

% ---------------------------------------------------------------------
% Run C version and plot Octave and C states and differences 
% ---------------------------------------------------------------------

system('../build_linux/unittest/tofdm');
load tofdm_out.txt;

stem_sig_and_error(1, 111, tx_bits_log_c, tx_bits_log - tx_bits_log_c, 'tx bits', [1 length(tx_bits_log) -1.5 1.5])
stem_sig_and_error(2, 211, real(tx_log_c), real(tx_log - tx_log_c), 'tx re', [1 length(tx_log_c) -0.1 0.1])
stem_sig_and_error(2, 212, imag(tx_log_c), real(tx_log - tx_log_c), 'tx im', [1 length(tx_log_c) -0.1 0.1])

% Run through checklist -----------------------------

check(W, W_c, 'W');
check(tx_bits_log, tx_bits_log_c, 'tx_bits');
check(tx_log, tx_log_c, 'tx');
