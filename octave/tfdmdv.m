% tfdmdv.m
%
% Octave script that tests the C port of the FDMDV modem.  This script loads
% the output of unittest/tfdmdv.c and compares it to the output of the
% reference versions of the same functions written in Octave.
%
% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

fdmdv; % load modem code
 
% Generate reference vectors using Octave implementation of FDMDV modem

global passes;
global fails;
passes = fails = 0;
frames = 25;
prev_tx_symbols = ones(Nc+1,1);
prev_rx_symbols = ones(Nc+1,1);
foff_phase_rect = 1;
track = 0;
fest_state = 0;

% Octave outputs we want to collect for comparison to C version

tx_bits_log = [];
tx_symbols_log = [];
tx_baseband_log = [];
tx_fdm_log = [];
pilot_baseband1_log = [];
pilot_baseband2_log = [];
pilot_lpf1_log = [];
pilot_lpf2_log = [];
S1_log = [];
S2_log = [];
foff_coarse_log = [];
foff_fine_log = [];
foff_log = [];
rx_baseband_log = [];
rx_filt_log = [];
env_log = [];
rx_timing_log = [];
rx_symbols_log = [];
rx_bits_log = []; 
sync_bit_log = [];  
track_log = [];

for f=1:frames

  % modulator

  tx_bits = get_test_bits(Nc*Nb);
  tx_bits_log = [tx_bits_log tx_bits];
  tx_symbols = bits_to_qpsk(prev_tx_symbols, tx_bits, 'dqpsk');
  prev_tx_symbols = tx_symbols;
  tx_symbols_log = [tx_symbols_log tx_symbols];
  tx_baseband = tx_filter(tx_symbols);
  tx_baseband_log = [tx_baseband_log tx_baseband];
  tx_fdm = fdm_upconvert(tx_baseband);
  tx_fdm_log = [tx_fdm_log tx_fdm];

  rx_fdm = real(tx_fdm);

  % demodulator

  [pilot prev_pilot pilot_lut_index prev_pilot_lut_index] = get_pilot(pilot_lut_index, prev_pilot_lut_index, M);

  [foff_coarse S1 S2] = rx_est_freq_offset(rx_fdm, pilot, prev_pilot, M);
  if track == 0
    foff = foff_coarse;
  end
  foff_log = [foff_log foff];
  foff_coarse_log = [foff_coarse_log foff_coarse];

  pilot_baseband1_log = [pilot_baseband1_log pilot_baseband1];
  pilot_baseband2_log = [pilot_baseband2_log pilot_baseband2];
  pilot_lpf1_log = [pilot_lpf1_log pilot_lpf1];
  pilot_lpf2_log = [pilot_lpf2_log pilot_lpf2];
  S1_log  = [S1_log S1];
  S2_log  = [S2_log S2];

  foff_rect = exp(j*2*pi*foff/Fs);

  for i=1:M
    foff_phase_rect *= foff_rect';
    rx_fdm_fcorr(i) = rx_fdm(i)*foff_phase_rect;
  end

  rx_baseband = fdm_downconvert(rx_fdm_fcorr, M);
  rx_baseband_log = [rx_baseband_log rx_baseband];

  rx_filt = rx_filter(rx_baseband, M);
  rx_filt_log = [rx_filt_log rx_filt];

  [rx_symbols rx_timing env] = rx_est_timing(rx_filt, rx_baseband, M);
  env_log = [env_log env];

  rx_timing_log = [rx_timing_log rx_timing];
  rx_symbols_log = [rx_symbols_log rx_symbols];

  [rx_bits sync_bit foff_fine] = qpsk_to_bits(prev_rx_symbols, rx_symbols, 'dqpsk');
  prev_rx_symbols = rx_symbols;
  rx_bits_log = [rx_bits_log rx_bits]; 
  foff_fine_log = [foff_fine_log foff_fine];
  sync_bit_log = [sync_bit_log sync_bit];  

  % freq est state machine

  [track fest_state] = freq_state(sync_bit, fest_state);
  track_log = [track_log track];
end

% Compare to the output from the C version

load ../unittest/tfdmdv_out.txt

% Helper functions to plot output of C verson and difference between Octave and C versions

function stem_sig_and_error(plotnum, subplotnum, sig, error, titlestr, axisvec)
  figure(plotnum)
  subplot(subplotnum)
  stem(sig);
  hold on;
  stem(error,'g');
  hold off;
  if nargin == 6
    axis(axisvec);
  end
  title(titlestr);
endfunction

function plot_sig_and_error(plotnum, subplotnum, sig, error, titlestr, axisvec)
  figure(plotnum)
  subplot(subplotnum)
  plot(sig);
  hold on;
  plot(error,'g');
  hold off;
  if nargin == 6
    axis(axisvec);
  end
  title(titlestr);
endfunction

% ---------------------------------------------------------------------------------------
% Plot output and test each C function
% ---------------------------------------------------------------------------------------

% fdmdv_get_test_bits() & bits_to_dqpsk_symbols()

n = 28;
stem_sig_and_error(1, 211, tx_bits_log_c(1:n), tx_bits_log(1:n) - tx_bits_log_c(1:n), 'tx bits', [1 n -1.5 1.5])
stem_sig_and_error(1, 212, real(tx_symbols_log_c(1:n/2)), real(tx_symbols_log(1:n/2) - tx_symbols_log_c(1:n/2)), 'tx symbols real', [1 n/2 -1.5 1.5])

% tx_filter()

diff = tx_baseband_log - tx_baseband_log_c;
c=3;
plot_sig_and_error(2, 211, real(tx_baseband_log_c(c,:)), real(sum(diff)), 'tx baseband real')
plot_sig_and_error(2, 212, imag(tx_baseband_log_c(c,:)), imag(sum(diff)), 'tx baseband imag')

% fdm_upconvert()

plot_sig_and_error(3, 211, real(tx_fdm_log_c), real(tx_fdm_log - tx_fdm_log_c), 'tx fdm real')
plot_sig_and_error(3, 212, imag(tx_fdm_log_c), imag(tx_fdm_log - tx_fdm_log_c), 'tx fdm imag')

% generate_pilot_lut()

plot_sig_and_error(4, 211, real(pilot_lut_c), real(pilot_lut - pilot_lut_c), 'pilot lut real')
plot_sig_and_error(4, 212, imag(pilot_lut_c), imag(pilot_lut - pilot_lut_c), 'pilot lut imag')

% rx_est_freq_offset()

plot_sig_and_error(5, 211, real(pilot_baseband1_log), real(pilot_baseband1_log - pilot_baseband1_log_c), 'pilot baseband1 real' )
plot_sig_and_error(5, 212, real(pilot_baseband2_log), real(pilot_baseband2_log - pilot_baseband2_log_c), 'pilot baseband2 real' )

plot_sig_and_error(6, 211, real(pilot_lpf1_log), real(pilot_lpf1_log - pilot_lpf1_log_c), 'pilot lpf1 real' )
plot_sig_and_error(6, 212, real(pilot_lpf2_log), real(pilot_lpf2_log - pilot_lpf2_log_c), 'pilot lpf2 real' )

plot_sig_and_error(7, 211, real(S1_log), real(S1_log - S1_log_c), 'S1 real' )
plot_sig_and_error(7, 212, imag(S1_log), imag(S1_log - S1_log_c), 'S1 imag' )

plot_sig_and_error(8, 211, real(S2_log), real(S2_log - S2_log_c), 'S2 real' )
plot_sig_and_error(8, 212, imag(S2_log), imag(S2_log - S2_log_c), 'S2 imag' )

plot_sig_and_error(9, 211, foff_coarse_log, foff_coarse_log - foff_coarse_log_c, 'Coarse Freq Offset' )
plot_sig_and_error(9, 212, foff_fine_log, foff_fine_log - foff_fine_log_c, 'Fine Freq Offset' )

plot_sig_and_error(10, 211, foff_log, foff_log - foff_log_c, 'Freq Offset' )
plot_sig_and_error(10, 212, track_log, track_log - track_log_c, 'Freq Track' )

c=15;
plot_sig_and_error(11, 211, real(rx_baseband_log(c,:)), real(rx_baseband_log(c,:) - rx_baseband_log_c(c,:)), 'Rx baseband real' )
plot_sig_and_error(11, 212, imag(rx_baseband_log(c,:)), imag(rx_baseband_log(c,:) - rx_baseband_log_c(c,:)), 'Rx baseband imag' )

plot_sig_and_error(12, 211, real(rx_filt_log(c,:)), real(rx_filt_log(c,:) - rx_filt_log_c(c,:)), 'Rx filt real' )
plot_sig_and_error(12, 212, imag(rx_filt_log(c,:)), imag(rx_filt_log(c,:) - rx_filt_log_c(c,:)), 'Rx filt imag' )

plot_sig_and_error(13, 211, env_log, env_log - env_log_c, 'env' )
plot_sig_and_error(13, 212, real(rx_symbols_log(c,:)), real(rx_symbols_log(c,:) - rx_symbols_log_c(c,:)), 'rx symbols' )

st=10*28;
en = 12*28;
plot_sig_and_error(14, 211, rx_timing_log, rx_timing_log - rx_timing_log_c, 'Rx Timing' )
stem_sig_and_error(14, 212, sync_bit_log_c, sync_bit_log - sync_bit_log_c, 'Sync bit', [1 n -1.5 1.5])

stem_sig_and_error(15, 211, rx_bits_log_c(st:en), rx_bits_log(st:en) - rx_bits_log_c(st:en), 'RX bits', [1 en-st -1.5 1.5])

% ---------------------------------------------------------------------------------------
% AUTOMATED CHECKS ------------------------------------------
% ---------------------------------------------------------------------------------------

function check(a, b, test_name)
  global passes;
  global fails;

  [m n] = size(a);
  printf("%s", test_name);
  for i=1:(25-length(test_name))
    printf(".");
  end
  printf(": ");  
  
  if abs(sum(a - b))/n < 1E-3
    printf("OK\n");
    passes++;
  else
    printf("FAIL\n");
    fails++;
  end
endfunction

check(tx_bits_log, tx_bits_log_c, 'tx_bits');
check(tx_symbols_log,  tx_symbols_log_c, 'tx_symbols');
check(tx_baseband_log, tx_baseband_log_c, 'tx_baseband');
check(tx_fdm_log, tx_fdm_log_c, 'tx_fdm');
check(pilot_lut, pilot_lut_c, 'pilot_lut');
check(pilot_baseband1_log, pilot_baseband1_log_c, 'pilot lpf1');
check(pilot_baseband2_log, pilot_baseband2_log_c, 'pilot lpf2');
check(S1_log, S1_log_c, 'S1');
check(S2_log, S2_log_c, 'S2');
check(foff_coarse_log, foff_coarse_log_c, 'foff_coarse');
check(foff_fine_log, foff_fine_log_c, 'foff_fine');
check(foff_log, foff_log_c, 'foff');
check(rx_baseband_log, rx_baseband_log_c, 'rx baseband');
check(rx_filt_log, rx_filt_log_c, 'rx filt');
check(env_log, env_log_c, 'env');
check(rx_timing_log, rx_timing_log_c, 'rx_timing');
check(rx_symbols_log, rx_symbols_log_c, 'rx_symbols');
check(rx_bits_log, rx_bits_log_c, 'rx bits');
check(sync_bit_log, sync_bit_log_c, 'sync bit');

printf("\npasses: %d fails: %d\n", passes, fails);
