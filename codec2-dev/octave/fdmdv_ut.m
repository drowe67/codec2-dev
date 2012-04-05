% fdmdv_ut.m
%
% Unit Test program for FDMDV modem.
%
% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

fdmdv;               % load modem code
 
% Simulation Parameters --------------------------------------

frames = 50;
EbNo_dB = 7.3;
Foff_hz = 0;
modulation = 'dqpsk';
hpa_clip = 150;

% ------------------------------------------------------------

tx_filt = zeros(Nc,M);
rx_symbols_log = [];
rx_phase_log = 0;
rx_timing_log = 0;
tx_pwr = 0;
noise_pwr = 0;
rx_fdm_log = [];
rx_baseband_log = [];
rx_bits_offset = zeros(Nc*Nb*2);
prev_tx_symbols = ones(Nc+1,1);
prev_rx_symbols = ones(Nc+1,1);
foff_log = [];
tx_baseband_log = [];
tx_fdm_log = [];

% BER stats

total_bit_errors = 0;
total_bits = 0;
bit_errors_log = [];
sync_log = [];
test_frame_sync_log = [];
test_frame_sync_state = 0;

% pilot states, used for copy of pilot at rx

pilot_rx_bit = 0;
pilot_symbol = sqrt(2);
pilot_freq = freq(Nc+1);
pilot_phase = 1;
pilot_filter_mem = zeros(1, Nfilter);
prev_pilot = zeros(M,1);

% fixed delay simuation

Ndelay = M+20;
rx_fdm_delay = zeros(Ndelay,1);

% ---------------------------------------------------------------------
% Eb/No calculations.  We need to work out Eb/No for each FDM carrier.
% Total power is sum of power in all FDM carriers
% ---------------------------------------------------------------------

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

% ---------------------------------------------------------------------
% Main loop 
% ---------------------------------------------------------------------

for f=1:frames

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
    % Time varying freq offset
    % Foff = Foff_hz + 100*sin(t*2*pi/(300*Fs));
    % t++;
    Foff = Foff_hz;
    freq_offset = exp(j*2*pi*Foff/Fs);
    phase_offset *= freq_offset;
    tx_fdm(i) = phase_offset*tx_fdm(i);
  end

  tx_fdm = real(tx_fdm);

  % HPA non-linearity

  tx_fdm(find(abs(tx_fdm) > hpa_clip)) = hpa_clip;
  tx_fdm_log = [tx_fdm_log tx_fdm];

  rx_fdm = tx_fdm;

  % AWGN noise

  noise = Ngain*randn(1,M);
  noise_pwr = 0.9*noise_pwr + 0.1*noise*noise'/M;
  rx_fdm += noise;
  rx_fdm_log = [rx_fdm_log rx_fdm];

  % Delay

  %rx_fdm_delay(1:Ndelay-M) = rx_fdm_delay(M+1:Ndelay);
  %rx_fdm_delay(Ndelay-M+1:Ndelay) = rx_fdm;
  rx_fdm_delay = rx_fdm;

  % -------------------
  % Demodulator
  % -------------------

  % frequency offset estimation and correction

  [pilot pilot_rx_bit pilot_symbol pilot_filter_mem pilot_phase] = generate_pilot_fdm(pilot_rx_bit, pilot_symbol, pilot_filter_mem, pilot_phase, pilot_freq);
  foff = rx_est_freq_offset(rx_fdm_delay, pilot, prev_pilot);
  prev_pilot = pilot;
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
  [rx_bits sync] = qpsk_to_bits(prev_rx_symbols, rx_symbols, modulation);
  prev_rx_symbols = rx_symbols;
  sync_log = [sync_log sync];

  % count bit errors if we find a test frame
  % Allow 15 frames for filter memories to fill and time est to settle

  [test_frame_sync bit_errors] = put_test_bits(rx_bits);
  if ((test_frame_sync == 1) && (f > 15))
    total_bit_errors = total_bit_errors + bit_errors;
    total_bits = total_bits + Ntest_bits;
    bit_errors_log = [bit_errors_log bit_errors];
  end
 
  % test frame sync state machine, just for more informative plots
    
  next_test_frame_sync_state = test_frame_sync_state;
  if (test_frame_sync_state == 0)
    if (test_frame_sync == 1)      
      next_test_frame_sync_state = 1;
      test_frame_count = 0;
    end
  end

  if (test_frame_sync_state == 1)
    % we only expect another test_frame_sync pulse every 4 symbols
    test_frame_count++;
    if (test_frame_count == 4)
      test_frame_count = 0;
      if ((test_frame_sync == 0))      
        next_test_frame_sync_state = 0;
      end
    end
  end
  test_frame_sync_state = next_test_frame_sync_state;
  test_frame_sync_log = [test_frame_sync_log test_frame_sync_state];
end

% ---------------------------------------------------------------------
% Print Stats
% ---------------------------------------------------------------------

ber = total_bit_errors / total_bits;

% Peak to Average Power Ratio from http://www.dsplog.com

papr = max(tx_fdm_log.*conj(tx_fdm_log)) / mean(tx_fdm_log.*conj(tx_fdm_log));
papr_dB = 10*log10(papr);

printf("Eb/No (meas): %2.2f (%2.2f) dB  %d bits  %d errors  BER: (%1.4f) PAPR: %1.2f dB  SNR: %2.1f dB\n", 
       EbNo_dB, 10*log10(0.25*tx_pwr*Fs/(Rs*Nc*noise_pwr)),
       total_bits, total_bit_errors, ber, papr_dB, SNR );

% ---------------------------------------------------------------------
% Plots
% ---------------------------------------------------------------------

figure(1)
clf;
[n m] = size(rx_symbols_log);
plot(real(rx_symbols_log(1:Nc+1,15:m)),imag(rx_symbols_log(1:Nc+1,15:m)),'+')
axis([-2 2 -2 2]);
title('Scatter Diagram');

figure(2)
clf;
subplot(211)
plot(rx_timing_log)
title('timing offset (samples)');
subplot(212)
plot(foff_log)
title('Freq offset (Hz)');

figure(3)
clf;
subplot(211)
plot(real(tx_fdm_log));
title('FDM Tx Signal');
subplot(212)
Nfft=Fs;
S=fft(rx_fdm_log,Nfft);
SdB=20*log10(abs(S));
plot(SdB(1:Fs/4))
title('FDM Rx Spectrum');

figure(4)
clf;
subplot(311)
stem(sync_log)
axis([0 frames 0 1.5]);
title('BPSK Sync')
subplot(312)
stem(bit_errors_log);
title('Bit Errors for test frames')
subplot(313)
plot(test_frame_sync_log);
axis([0 frames 0 1.5]);
title('Test Frame Sync')

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

