% fdmdv_ut_freq_off.m
% David Rowe 17 June 2014
%
% Unit Test program for freq offset estimation in FDMDV modem.  
%
% Copyright David Rowe 2012 This program is
% distributed under the terms of the GNU General Public License
% Version 2

% [ ] sweep of different delays
% [ ] sweep of Eb/No
% [ ] sweep of freq offsets
% [ ] plot/print pass fail/relevant stats
%      + variance
%      + histogram of freq ests?

fdmdv;  % load modem code
 
% Simulation Parameters --------------------------------------

frames = 100;
EbNovec = 10;
Foff_hz = 0;

% ---------------------------------------------------------------------
% Eb/No calculations.  We need to work out Eb/No for each FDM carrier.
% Total power is sum of power in all FDM carriers
% ---------------------------------------------------------------------

function Nsd = calc_Nsd_from_EbNo(EbNo_dB)
  global Rs;
  global Nb;
  global Fs;

  C = 1; % power of each FDM carrier (energy/sample).  Total Carrier power should = Nc*C = Nc
  N = 1; % total noise power (energy/sample) of noise source across entire bandwidth

  % Eb  = Carrier power * symbol time / (bits/symbol)
  %     = C *(1/Rs) / Nb
  Eb_dB = 10*log10(C) - 10*log10(Rs) - 10*log10(Nb);

  No_dBHz = Eb_dB - EbNo_dB;

  % Noise power = Noise spectral density * bandwidth
  % Noise power = Noise spectral density * Fs/2 for real signals
  N_dB = No_dBHz + 10*log10(Fs/2);
  Ngain_dB = N_dB - 10*log10(N);
  Nsd = 10^(Ngain_dB/20);
end

% ------------------------------------------------------------

modulation = 'dqpsk';
tx_filt = zeros(Nc,M);
tx_fdm_log = [];
rx_fdm_log = [];
prev_tx_symbols = ones(Nc+1,1);
ferr = 0;
foff = 0;
foff_log = [];

% fixed delay simuation

Ndelay = M+20;
rx_fdm_delay = zeros(Ndelay,1);

% ---------------------------------------------------------------------
% Eb/No calculations.  We need to work out Eb/No for each FDM carrier.
% Total power is sum of power in all FDM carriers
% ---------------------------------------------------------------------

% freq offset simulation states

phase_offset = 1;
freq_offset = exp(j*2*pi*Foff_hz/Fs);
foff_phase = 1;
t = 0;
foff = 0;

% ---------------------------------------------------------------------
% Main loop 
% ---------------------------------------------------------------------

for ne = 1:length(EbNovec)
   EbNo_dB = EbNovec(ne);
   Nsd = calc_Nsd_from_EbNo(EbNo_dB);

  for f=1:frames

    % ------------------- Modulator -------------------

    tx_bits = get_test_bits(Nc*Nb); 
    tx_symbols = bits_to_psk(prev_tx_symbols, tx_bits, modulation); 
    prev_tx_symbols = tx_symbols; 
    tx_baseband = tx_filter(tx_symbols); 
    tx_fdm = fdm_upconvert(tx_baseband);
    tx_fdm_log = [tx_fdm_log real(tx_fdm)];

    % ------------------- Channel simulation -------------------

    % frequency offset

    Foff = Foff_hz; for i=1:M freq_offset = exp(j*2*pi*Foff/Fs);
    phase_offset *= freq_offset; tx_fdm(i) = phase_offset*tx_fdm(i); end

     rx_fdm = real(tx_fdm);

    % AWGN noise

    noise = Nsd*randn(1,M); 
    rx_fdm += noise; 
    rx_fdm_log = [rx_fdm_log rx_fdm];

    % Delay

    rx_fdm_delay(1:Ndelay-M) = rx_fdm_delay(M+1:Ndelay);
    rx_fdm_delay(Ndelay-M+1:Ndelay) = rx_fdm; 
    %rx_fdm_delay = rx_fdm;
  
    % ------------------- Freq Offset Est -------------------

    % frequency offset estimation and correction, need to call
    % rx_est_freq_offset even in track mode to keep states updated

    [pilot prev_pilot pilot_lut_index prev_pilot_lut_index] = ...
    get_pilot(pilot_lut_index, prev_pilot_lut_index, M); 
    [foff_coarse S1 S2] = rx_est_freq_offset(rx_fdm_delay, pilot, prev_pilot, M);
    foff_log(ne,f) = foff_coarse;
  end
end

% ---------------------------------------------------------------------
% Print Stats
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% Plots
% ---------------------------------------------------------------------

figure(1)
clf
for ne = 1:length(EbNovec)
  foff_std(ne) = std(foff_log(ne,:));
end
plot(EbNovec,foff_std)
xlabel("Eb/No (dB)")
ylabel("Std Dev")

figure(2)
clf;
plot(foff_log(1,:))
xlabel("Frames")
ylabel("Freq offset estimate")

figure(3)
clf;
hist(foff_log(1,:))
