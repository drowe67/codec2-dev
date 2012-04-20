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
tx_bits_log = [];
tx_symbols_log = [];
tx_baseband_log = [];
tx_fdm_log = [];

for f=1:frames
  tx_bits = get_test_bits(Nc*Nb);
  tx_bits_log = [tx_bits_log tx_bits];
  tx_symbols = bits_to_qpsk(prev_tx_symbols, tx_bits, 'dqpsk');
  prev_tx_symbols = tx_symbols;
  tx_symbols_log = [tx_symbols_log tx_symbols];
  tx_baseband = tx_filter(tx_symbols);
  tx_baseband_log = [tx_baseband_log tx_baseband];
  tx_fdm = fdm_upconvert(tx_baseband);
  tx_fdm_log = [tx_fdm_log tx_fdm];
end

% Compare to the output from the C version

load ../unittest/tfdmdv_out.txt

figure(1)
subplot(211)
n = 28;
stem(tx_bits_log_c(1:n));
hold on;
stem(tx_bits_log(1:n) - tx_bits_log_c(1:n),'g');
hold off;
axis([1 n -1.5 1.5])
title('tx bits')
subplot(212)
stem(real(tx_symbols_log_c(1:n/2)));
hold on;
stem(tx_symbols_log(1:n/2) - tx_symbols_log_c(1:n/2),'g');
hold off;
axis([1 n/2 -1.5 1.5])
title('tx symbols real')

figure(2)
clf;
diff = tx_baseband_log - tx_baseband_log_c;
subplot(211)
c=3;
plot(real(tx_baseband_log_c(c,:)));
hold on;
plot(real(sum(diff)),'g')
hold off;
title('tx baseband real')
subplot(212)
plot(imag(tx_baseband_log_c(c,:)));
hold on;
plot(imag(sum(diff)),'g')
hold off;
title('tx baseband imag')

figure(3)
clf
subplot(211)
plot(real(tx_fdm_log_c));
hold on;
plot(real(tx_fdm_log - tx_fdm_log_c),'g');
hold off;
title('tx fdm real')
subplot(212)
plot(imag(tx_fdm_log_c));
hold on;
plot(imag(tx_fdm_log - tx_fdm_log_c),'g');
hold off;
title('tx fdm imag')

figure(4)
clf
subplot(211)
plot(real(pilot_lut_c));
hold on;
plot(real(pilot_lut - pilot_lut_c),'g');
hold off;
title('pilot lut real')
subplot(212)
plot(imag(pilot_lut_c));
hold on;
plot(imag(pilot_lut - pilot_lut_c),'g');
hold off;
title('pilot lut imag')

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

printf("\npasses: %d fails: %d\n", passes, fails);
