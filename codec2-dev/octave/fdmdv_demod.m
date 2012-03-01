% fdmdv_demod.m
%
% Demodulator function for FDMDV modem.
%
% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

function fdmdv_demod(rawfilename)

fdmdv; % include modem code

fin = fopen(rawfilename, "rb");
rx_fdm = fread(fin, Inf, "short");
gain = 1000;
rx_fdm /= gain;
frames = floor(length(rx_fdm)/M);

total_bit_errors = 0;
total_bits = 0;

rx_timing_log = [];
rx_symbols_log = [];

% Main loop ----------------------------------------------------

for i=1:frames
  rx_baseband = fdm_downconvert(rx_fdm((i-1)*M+1:i*M));
  rx_filt = rx_filter(rx_baseband);

  [rx_symbols rx_timing] = rx_est_timing(rx_filt, rx_baseband);
  rx_timing_log = [rx_timing_log rx_timing];

  rx_symbols_log = [rx_symbols_log rx_symbols];
  rx_bits = qpsk_to_bits(rx_symbols);

  [sync bit_errors] = put_test_bits(rx_bits);

  if sync == 1
    total_bit_errors = total_bit_errors + bit_errors;
    total_bits = total_bits + Ntest_bits;
  endif

  end

end

ber = total_bit_errors/total_bits;
printf("%d bits  %d errors  Meas BER: %1.4f\n", total_bits, total_bit_errors,ber);

figure(1)
clf;
[n m] = size(rx_symbols_log);
plot(real(rx_symbols_log(:,20:m)),imag(rx_symbols_log(:,20:m)),'+')
figure(2)
clf;
subplot(211)
plot(rx_timing_log)
subplot(212)
Nfft=Fs;
S=fft(rx_fdm,Nfft);
SdB=20*log10(abs(S));
plot(SdB(1:Fs/4))
