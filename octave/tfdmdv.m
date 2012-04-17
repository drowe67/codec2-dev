% tfdmdv.m
%
% Octave script that evaluates the output of tfdmdv.c Unit Test program for FDMDV modem.
%
% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

fdmdv; % load modem code
 
% Generate reference vectors using Octave implementation of FDMDV modem

passes = fails = 0;
frames = 10;
prev_tx_symbols = ones(Nc+1,1);
tx_bits_log = [];
tx_symbols_log = [];

for f=1:frames

  tx_bits = get_test_bits(Nc*Nb);
  tx_bits_log = [tx_bits_log tx_bits];
  tx_symbols = bits_to_qpsk(prev_tx_symbols, tx_bits, 'dqpsk');
  prev_tx_symbols = tx_symbols;
  tx_symbols_log = [tx_symbols_log tx_symbols];

end

% Compare to the output from the C version

load ../unittest/tfdmdv_out.txt

figure(1)
subplot(211)
plot(tx_bits_log - tx_bits_tfdmdv);
title('tx bits')
subplot(212)
plot(tx_symbols_log - tx_symbols_tfdmdv);
title('tx symbols')

if sum(tx_bits_log - tx_bits_tfdmdv) == 0
  printf("fdmdv_get_test_bits..: OK\n");
  passes++;
else;
  printf("fdmdv_get_test_bits..: FAIL\n");
  fails++;
end
 
if sum(tx_symbols_log - tx_symbols_tfdmdv) == 0
  printf("bits_to_dqpsk_symbols: OK\n");
  passes++;
else;
  printf("bits_to_dqpsk_symbols: FAIL\n");
  fails++;
end

printf("\npasses: %d fails: %d\n", passes, fails);
