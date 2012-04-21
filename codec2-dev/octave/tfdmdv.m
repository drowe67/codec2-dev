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

passes = fails = 0;
frames = 25;
prev_tx_symbols = ones(Nc+1,1);

% Octave outputs we want to collect for comparison to C version

tx_bits_log = [];
tx_symbols_log = [];
tx_baseband_log = [];
tx_fdm_log = [];
pilot_baseband1_log = [];
pilot_baseband2_log = [];
pilot_lpf1_log = [];
pilot_lpf2_log = [];
s1_log = [];
s2_log = [];

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

  [foff s1 s2] = rx_est_freq_offset(rx_fdm, pilot, prev_pilot, M);

  pilot_baseband1_log = [pilot_baseband1_log pilot_baseband1];
  pilot_baseband2_log = [pilot_baseband2_log pilot_baseband2];
  pilot_lpf1_log = [pilot_lpf1_log pilot_lpf1];
  pilot_lpf2_log = [pilot_lpf2_log pilot_lpf2];
  s1_log  = [s1_log s1];
  s2_log  = [s2_log s2];

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
stem_sig_and_error(1, 212, real(tx_symbols_log_c(1:n/2)), tx_symbols_log(1:n/2) - tx_symbols_log_c(1:n/2), 'tx symbols real', [1 n/2 -1.5 1.5])

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

plot_sig_and_error(7, 211, real(s1_log), real(s1_log - s1_log_c), 's1 real' )
plot_sig_and_error(7, 212, real(s2_log), real(s2_log - s2_log_c), 's2 real' )

% ---------------------------------------------------------------------------------------
% AUTOMATED CHECKS ------------------------------------------
% ---------------------------------------------------------------------------------------

if sum(tx_bits_log - tx_bits_log_c) == 0
  printf("fdmdv_get_test_bits..: OK\n");
  passes++;
else;
  printf("fdmdv_get_test_bits..: FAIL\n");
  fails++;
end
 
if sum(tx_symbols_log - tx_symbols_log_c) == 0
  printf("bits_to_dqpsk_symbols: OK\n");
  passes++;
else;
  printf("bits_to_dqpsk_symbols: FAIL\n");
  fails++;
end

if sum(tx_baseband_log - tx_baseband_log_c) < 1E-3
  printf("tx_filter............: OK\n");
  passes++;
else;
  printf("tx_filter............: FAIL\n");
  fails++;
end

if sum(tx_fdm_log - tx_fdm_log_c) < 1E-3
  printf("tx_fdm...............: OK\n");
  passes++;
else;
  printf("tx_fdm...............: FAIL\n");
  fails++;
end

if sum(pilot_lut - pilot_lut_c) < 1E-3
  printf("pilot_lut............: OK\n");
  passes++;
else;
  printf("pilot_lut............: FAIL\n");
  fails++;
end

[m n] = size(pilot_baseband1_log);
if sum(pilot_baseband1_log - pilot_baseband1_log_c)/n < 1E-3
  printf("pilot_baseband1_log..: OK\n");
  passes++;
else;
  printf("pilot_baseband1_log..: FAIL\n");
  fails++;
end

[m n] = size(pilot_baseband2_log);
if sum(pilot_baseband2_log - pilot_baseband2_log_c)/n < 1E-3
  printf("pilot_baseband2_log..: OK\n");
  passes++;
else;
  printf("pilot_baseband2_log..: FAIL\n");
  fails++;
end

printf("\npasses: %d fails: %d\n", passes, fails);
