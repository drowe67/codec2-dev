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
 
% ---------------------------------------------------------------------
% Eb/No calculations.  We need to work out Eb/No for each FDM carrier.
% Total power is sum of power in all FDM carriers
% ---------------------------------------------------------------------

function [Nsd SNR] = calc_Nsd_from_EbNo(EbNo_dB)
  global Rs;
  global Nb;
  global Nc;
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

  % C/No = Carrier Power/noise spectral density
  %      = power per carrier*number of carriers / noise spectral density
  CNo_dB = 10*log10(C)  + 10*log10(Nc) - No_dBHz;

  % SNR in equivalent 3000 Hz SSB channel

  B = 3000;
  SNR = CNo_dB - 10*log10(B);
end

% ------------------------------------------------------------

function sim_out = freq_off_est_test(sim_in)
  global Nc;
  global Nb;
  global M;
  global Fs;
  global pilot_lut_index;
  global prev_pilot_lut_index;
  global pilot_lpf1;
  global Npilotlpf;

  EbNovec = sim_in.EbNovec;
  Ndelay = sim_in.delay;
  frames = sim_in.frames;
  startup_delay = sim_in.startup_delay;
  allowable_error = sim_in.allowable_error;
  foff_hz = sim_in.foff_hz;

  % ---------------------------------------------------------------------
  % Main loop 
  % ---------------------------------------------------------------------

  for ne = 1:length(EbNovec)
     EbNo_dB = EbNovec(ne);
     [Nsd SNR] = calc_Nsd_from_EbNo(EbNo_dB);
     hits = 0;

     tx_filt = zeros(Nc,M);
     prev_tx_symbols = ones(Nc+1,1);

     tx_fdm_log = [];
     rx_fdm_log = [];
     pilot_lpf1_log = [];
     S1_log = [];
     rx_fdm_delay = zeros(M+Ndelay,1);

     % freq offset simulation states

     phase_offset = 1;
     freq_offset = exp(j*2*pi*foff_hz/Fs);
     foff_phase = 1;

     for f=1:frames

      % ------------------- Modulator -------------------

      tx_bits = get_test_bits(Nc*Nb); 
      tx_symbols = bits_to_psk(prev_tx_symbols, tx_bits, 'dqpsk'); 
      prev_tx_symbols = tx_symbols; 
      tx_baseband = tx_filter(tx_symbols); 
      tx_fdm = fdm_upconvert(tx_baseband);
      tx_fdm_log = [tx_fdm_log real(tx_fdm)];

      % ------------------- Channel simulation -------------------

      % frequency offset

      foff = foff_hz; 
      for i=1:M 
        freq_offset = exp(j*2*pi*foff/Fs);
        phase_offset *= freq_offset; 
        tx_fdm(i) = phase_offset*tx_fdm(i); 
      end

      rx_fdm = real(tx_fdm);

      % AWGN noise

      noise = Nsd*randn(1,M); 
      rx_fdm += noise; 
      rx_fdm_log = [rx_fdm_log rx_fdm];

      % Delay

      rx_fdm_delay(1:Ndelay) = rx_fdm_delay(M+1:M+Ndelay);
      rx_fdm_delay(Ndelay+1:M+Ndelay) = rx_fdm; 

      % ------------------- Freq Offset Est -------------------

      % frequency offset estimation and correction, need to call
      % rx_est_freq_offset even in track mode to keep states updated

      [pilot prev_pilot pilot_lut_index prev_pilot_lut_index] = ...
          get_pilot(pilot_lut_index, prev_pilot_lut_index, M); 
      [foff_coarse S1 S2] = rx_est_freq_offset(rx_fdm_delay, pilot, prev_pilot, M);
      pilot_lpf1_log = [pilot_lpf1_log pilot_lpf1(Npilotlpf-M+1:Npilotlpf)];
      S1_log(f,:) = fftshift(S1);

      foff_log(ne,f) = foff_coarse;

      if (f > startup_delay) && (abs(foff_coarse < foff_hz) < allowable_error)
        hits++;
      end
    end

    % results for this EbNo value

    sim_out.foff_sd(ne) = std(foff_log(ne,startup_delay:frames));
    sim_out.hits = hits;
    sim_out.hits_percent = 100*sim_out.hits/(frames-startup_delay);

    printf("EbNo (dB): %3.2f  SNR (3kHz dB): %3.2f  std dev (Hz): %3.2f  Hits: %d (%3.2f%%)\n", ...
           EbNo_dB, SNR, sim_out.foff_sd(ne), sim_out.hits, sim_out.hits_percent);

    % plots if single dimension vector

    if length(EbNovec) == 1
      figure(2)
      clf;
      plot(foff_log(ne,:))
      xlabel("Frames")
      ylabel("Freq offset estimate")

      figure(3)
      clf;
      hist(foff_log(ne,:));

      figure(4)
      [n m] = size(S1_log);
      mesh(-200+400*(0:m-1)/256,1:n,abs(S1_log(:,:)))
    end
  end
end

% ---------------------------------------------------------------------
% Run Automated Tests
% ---------------------------------------------------------------------

sim_in.EbNovec = 0:10;
sim_in.delay = M/2;
sim_in.frames = 20;
sim_in.foff_hz = 0;
sim_in.startup_delay = 10;
sim_in.allowable_error = 5;

sim_out = freq_off_est_test(sim_in);

figure(1)
clf
plot(sim_in.EbNovec,sim_out.foff_sd)
xlabel("Eb/No (dB)")
ylabel("Std Dev")

