% fdmdv_ut.m
%
% Unit Test program for FDMDV modem.
%
% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

clear all;

fdmdv;               % load modem code
 
% Simulation Parameters --------------------------------------

frames = 50;
EbNo_dB = 7.3;
Foff_hz = 0;
modulation = 'dqpsk';

% ------------------------------------------------------------

tx_filt = zeros(Nc,M);
rx_symbols_log = zeros(Nc,1);
rx_phase_log = 0;
rx_timing_log = 0;
tx_pwr = 0;
noise_pwr = 0;
total_bit_errors = 0;
total_bits = 0;
rx_fdm_log = [];
rx_baseband_log = [];
rx_bits_offset = zeros(Nc*Nb*2);
prev_tx_symbols = sqrt(2)*ones(Nc,1)*exp(j*pi/4);
prev_rx_symbols = sqrt(2)*ones(Nc,1)*exp(j*pi/4);
foff_log = [];
tx_baseband_log = [];

Ndelay = M+20;
rx_fdm_delay = zeros(Ndelay,1);

% Eb/No calculations.  We need to work out Eb/No for each FDM carrier.
% Total power is sum of power in all FDM carriers

C = 1; % power of each FDM carrier (energy/sample).  Total Carrier power should = Nc*C = Nc
N = 1; % total noise power (energy/sample) of noise source across entire bandwidth

% Eb  = Carrier power * symbol time / (bits/symbol)
%     = C *(1/Rs) / 2
Eb_dB = 10*log10(C) - 10*log10(Rs) - 10*log10(2);

No_dBHz = Eb_dB - EbNo_dB;

% Noise power = Noise spectral density * bandwidth
% Noise power = Noise spectral density * Fs/2 for real signals
N_dB = No_dBHz + 10*log10(Fs/2);
Ngain_dB = N_dB - 10*log10(N);
Ngain = 10^(Ngain_dB/20);

% C/No = Carrier Power/noise spectral density
%      = power per carrier*number of carriers / noise spectral density
CNo_dB = 10*log10(C)  + 10*log10(Nc) - No_dBHz;

% SNR in equivalent 3000 Hz SSB channel

B = 3000;
SNR = CNo_dB - 10*log10(B);

phase_offset = 1;
freq_offset = exp(j*2*pi*Foff_hz/Fs);
foff_phase = 1;
t = 0;

% Main loop ----------------------------------------------------

for i=1:frames

  % -------------------
  % Modulator
  % -------------------

  tx_bits = get_test_bits(Nc*Nb);
  tx_symbols = bits_to_qpsk(prev_tx_symbols, tx_bits, modulation);
  prev_tx_symbols = tx_symbols;
  tx_baseband = tx_filter(tx_symbols);
  tx_baseband_log = [tx_baseband_log tx_baseband];
  tx_fdm = fdm_upconvert(tx_baseband);
  tx_pwr = 0.9*tx_pwr + 0.1*real(tx_fdm)*real(tx_fdm)'/(M);

  % -------------------
  % Channel simulation
  % -------------------

  % frequency offset

  for i=1:M
    Foff = Foff_hz + 100*sin(t*2*pi/(300*Fs));
    t++;
    freq_offset = exp(j*2*pi*Foff/Fs);
    phase_offset *= freq_offset;
    rx_fdm(i) = phase_offset*real(tx_fdm(i));
  end

  % AWGN noise

  noise = Ngain*randn(1,M);
  noise_pwr = 0.9*noise_pwr + 0.1*noise*noise'/M;
  rx_fdm += noise;
  rx_fdm_log = [rx_fdm_log rx_fdm];

  % delay

  rx_fdm_delay(1:Ndelay-M) = rx_fdm_delay(M+1:Ndelay);
  rx_fdm_delay(Ndelay-M+1:Ndelay) = rx_fdm;

  % -------------------
  % Demodulator
  % -------------------

  % frequency offset estimation and correction

  foff = rx_est_freq_offset(rx_fdm);
  foff_log = [ foff_log foff ];
  %foff = 0;
  foff_rect = exp(j*2*pi*foff/Fs);

  for i=1:M
    foff_phase *= foff_rect';
    rx_fdm_delay(i) = rx_fdm_delay(i)*foff_phase;
  end

  % baseband processing

  rx_baseband = fdm_downconvert(rx_fdm_delay(1:M));
  rx_baseband_log = [rx_baseband_log rx_baseband];
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

  % count bit errors

  [sync bit_errors] = put_test_bits(rx_bits);
  if (sync == 1)
    total_bit_errors = total_bit_errors + bit_errors;
    total_bits = total_bits + Ntest_bits;
  end

end

ber = total_bit_errors/total_bits;
printf("Eb/No (meas): %2.2f (%2.2f) dB  %d bits  %d errors  QPSK BER (meas): %1.4f (%1.4f)\n", 
       EbNo_dB, 10*log10(0.25*tx_pwr*Fs/(Rs*Nc*noise_pwr)),
       total_bits, total_bit_errors, 0.5*erfc(sqrt(10.^(EbNo_dB/10))), ber );

figure(1)
clf;
[n m] = size(rx_symbols_log);
plot(real(rx_symbols_log(:,20:m)),imag(rx_symbols_log(:,20:m)),'+')

figure(2)
clf;
subplot(211)
plot(rx_timing_log)
title('timing offset (samples)');
subplot(212)
plot(foff_log)
title('Freq offset (Hz)');

%figure(3)
%clf;
%Nfft=Fs;
%S=fft(rx_fdm_log,Nfft);
%SdB=20*log10(abs(S));
%plot(-Fs/2+1:Fs/2,fftshift(SdB))
%plot(SdB(1:Fs/4))



% TODO
%   + handling sample slips, extra plus/minus samples
%   + simulating sample clock offsets
%   + timing, frequency offset get function
%   + memory of recent tx and rx signal for spectrogram
%   + scatter diagram get function

% DPSK
% add phase offset, frequency offset, timing offset
% fading simulator
% file I/O to test with Codec
% code to measure sync recovery time
% dump file type plotting & instrumentation
% determine if error pattern is bursty
% HF channel simulation
% Offset or pi/4 QPSK and tests with real tx HPA
% real time SNR get function
%
% phase estimator not working too well and would need a UW
% to resolve ambiguity.  But this is probably worth it for
% 3dB.  Test with small freq offset

% Implementation loss BER issues:
%   QPSK mapping
%   Single sided noise issues
%   interference between carriers due to excess BW
%   Crappy RN coeffs
%   timing recovery off by one
%   Use a true PR sequence
%   Sensitivity to Fs
%   Try BPSK
%   second term of QPSK eqn

% Faster sync:
%
% 1/ Maybe Ask Bill how we can sync in less than 10 symbols?  What to
% put in filter memories?  If high SNR should sync very quickly
% Maybe start feeding in symbols to middle of filter memory?  Or do timing
% sync before rx filtering at start?

