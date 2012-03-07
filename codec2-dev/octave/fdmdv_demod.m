% fdmdv_demod.m
%
% Demodulator function for FDMDV modem.
%
% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

function fdmdv_demod(rawfilename, nbits)

fdmdv; % include modem code

modulation = 'dqpsk';
Foff_hz = 0;
phase_offset = 1;
freq_offset = exp(j*2*pi*Foff_hz/Fs);

fin = fopen(rawfilename, "rb");
rx_fdm = fread(fin, Inf, "short");
gain = 1000;
rx_fdm /= gain;
if (nargin == 1)
  frames = floor(length(rx_fdm)/M);
else
    frames = nbits/(Nc*Nb);
endif

prev_rx_symbols = sqrt(2)*ones(Nc,1)*exp(j*pi/4);

total_bit_errors = 0;
total_bits = 0;

rx_timing_log = [];
rx_symbols_log = [];
rx_phase_log = [];
sync_log = [];

% Main loop ----------------------------------------------------

for i=1:frames
  for n=1:M
    phase_offset *= freq_offset;
    rx_fdm((i-1)*M+1+n) *= phase_offset;
  end
  rx_baseband = fdm_downconvert(rx_fdm((i-1)*M+1:i*M));
  rx_filt = rx_filter(rx_baseband);

  [rx_symbols rx_timing] = rx_est_timing(rx_filt, rx_baseband);
  rx_timing_log = [rx_timing_log rx_timing];

  %rx_phase = rx_est_phase(rx_symbols);
  %rx_phase_log = [rx_phase_log rx_phase];
  %rx_symbols = rx_symbols*exp(j*rx_phase);


  if strcmp(modulation,'dqpsk')
    rx_symbols_log = [rx_symbols_log rx_symbols.*conj(prev_rx_symbols)*exp(j*pi/4)];
  else
    rx_symbols_log = [rx_symbols_log rx_symbols];
  endif
  rx_bits = qpsk_to_bits(prev_rx_symbols, rx_symbols, modulation);
  prev_rx_symbols = rx_symbols;

  [sync bit_errors] = put_test_bits(rx_bits);
  sync_log = [sync_log sync];

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
%plot(rx_fdm);
figure(2)
clf;
%subplot(211)
plot(rx_timing_log)
%subplot(212)
%plot(sync_log)
figure(3)
clf;
for i=1:Nc
  subplot(4,4,i);
  mx = max(abs(rx_symbols_log(i,20:m)));
  plot(real(rx_symbols_log(i,20:m)),imag(rx_symbols_log(i,20:m)),'+')
  axis([-mx mx -mx mx]);
end
