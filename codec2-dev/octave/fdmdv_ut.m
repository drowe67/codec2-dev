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
rand('state',1); 
randn('state',1);
 
frames = 50;
tx_filt = zeros(Nc,M);
rx_symbols_log = zeros(Nc,1);
rx_phase_log = 0;
rx_timing_log = 0;
tx_pwr = 0;
noise_pwr = 0;
total_bit_errors = 0;
total_bits = 0;
rx_fdm_log = [];

% Eb/No calculations.  We need to work out Eb/No for each FDM carrier.
% Total power is sum of power in all FDM carriers

C = 1;  % power of each FDM carrier (energy/sample)
N = 1;  % total noise power (energy/sample) of noise source before scaling 
        % by Ngain

EbNo_dB = 40;

% Eb  = Carrier power * symbol time / (bits/symbol)
%     = C *(Rs/Fs) / 2
Eb_dB = 10*log10(C) + 10*log10(Rs) - 10*log10(Fs) - 10*log10(2);

No_dBHz = Eb_dB - EbNo_dB;

% Noise power = Noise spectral density * bandwidth;
N_dB = No_dBHz + 10*log10(Fs);
Ngain_dB = N_dB - 10*log10(N);
Ngain = 10^(Ngain_dB/20);

% C/No = Carrier Power/noise spectral denity
%      = power per carrier*number of carriers / noise spectral denity
CNo_dB = 10*log10(C)  + 10*log10(Nc) - No_dBHz;

% SNR in equivalent 2400 Hz SSB channel

B = 2400;
SNR = CNo_dB - 10*log10(B);


% Main loop ----------------------------------------------------

for i=1:frames
  tx_bits = get_test_bits(Nc*Nb);
  tx_symbols = bits_to_qpsk(tx_bits);
  tx_baseband = tx_filter(tx_symbols);
  tx_fdm = fdm_upconvert(tx_baseband);
  tx_pwr = 0.9*tx_pwr + 0.1*tx_fdm*tx_fdm'/(M);

  noise = Ngain/sqrt(2)*[randn(1,M) + j*randn(1,M)];
  noise_pwr = 0.9*noise_pwr + 0.1*noise*noise'/M;
  rx_fdm = tx_fdm + noise;
  rx_fdm_log = [rx_fdm_log rx_fdm];

  rx_baseband = fdm_downconvert(rx_fdm);
  rx_filt = rx_filter(rx_baseband);

  [rx_symbols rx_timing] = rx_est_timing(rx_filt, rx_baseband);
  rx_timing_log = [rx_timing_log rx_timing];

  %rx_phase = rx_est_phase(rx_symbols);
  %rx_phase_log = [rx_phase_log rx_phase];
  %rx_symbols = rx_symbols*exp(j*rx_phase);

  rx_symbols_log = [rx_symbols_log rx_symbols];
  rx_bits = qpsk_to_bits(rx_symbols);

  [sync bit_errors] = put_test_bits(rx_bits);
  if (sync == 1)
    total_bit_errors = total_bit_errors + bit_errors;
    total_bits = total_bits + Ntest_bits;
  end

end

ber = total_bit_errors/total_bits;
printf("Eb/No: %2.2f dB  %d bits  %d errors  Meas BER: %1.4f  Theor BER: %1.4f\n", EbNo_dB, 
      total_bits, total_bit_errors, ber, 0.5*erfc(sqrt(10.^(EbNo_dB/10))));

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
S=fft(rx_fdm_log,Nfft);
SdB=20*log10(abs(S));
%plot(-Fs/2+1:Fs/2,fftshift(SdB))
plot(SdB(1:Fs/4))

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
