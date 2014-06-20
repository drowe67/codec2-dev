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
% [ ] step change in foff
%     + time to respond
% [ ] plot/print pass fail/relevant stats
%      + variance
%      + histogram of freq ests?

fdmdv;  % load modem code
hf_sim; % load hf sim code

% ---------------------------------------------------------------------
% Eb/No calculations.  We need to work out Eb/No for each FDM carrier.
% Total power is sum of power in all FDM carriers.  These calcs set the
% Eb/No of the data carriers, Eb/No of pilot will be higher.
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
  CNo_dB = 10*log10(C) + 10*log10(Nc) - No_dBHz;

  % SNR in equivalent 3000 Hz SSB channel, adding extra power for pilot to get
  % true SNR.

  B = 3000;
  SNR = CNo_dB - 10*log10(B) + 10*log10((Nc+4)/Nc);
end

% we keep a m sample buffer in sample_memory
% update sample_memory with n samples each time this function is called
% outputs one nfft2 slice of spectrogram in dB.  Good idea to make nfft2 a power of 2

function [S, states_out] = spectrogram_update(samples, n, states_in)
  sample_memory = states_in.sample_memory;
  m             = states_in.m;
  nfft2         = states_in.nfft2;
  lower_clip_dB = states_in.lower_clip_dB;
  dec           = states_in.dec;

  sample_memory(1:m-n)   = sample_memory(n+1:m);
  sample_memory(m-n+1:m) = samples;

  F = fft(sample_memory .* hanning(m)', 2*nfft2);
  S = 20*log10(abs(F(1:dec:nfft2))/(nfft2));
  S(find(S < lower_clip_dB)) = lower_clip_dB;    % clip lower limit

  states_out = states_in;
  states_out.sample_memory = sample_memory;
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
  global spread;
  global spread_2ms;
  global hf_gain;

  EbNovec  = sim_in.EbNovec;
  Ndelay   = sim_in.delay;
  frames   = sim_in.frames;
  startup_delay = sim_in.startup_delay;
  allowable_error = sim_in.allowable_error;
  foff_hz  = sim_in.foff_hz;
  hf_sim   = sim_in.hf_sim;
  hf_delay = floor(sim_in.hf_delay_ms*Fs/1000);

  % work out gain for HF model
  % e = sum((g*s)^2) = g*g*sum(s^2) = N, g = sqrt(N/sum(s^2))
  % compute so e=N

  s1 = spread(1:frames*M);
  s2 = [zeros(hf_delay,1); spread_2ms(1:frames*M)];
  s2 = s2(1:frames*M);

  p = (s1+s2)'*(s1+s2);
  hf_gain = sqrt(frames*M/p);
  p2 = (hf_gain*(s1+s2))'*(hf_gain*(s1+s2));

  % spectrogram states

  spec_states.m             = 4*M;
  spec_states.nfft2         = 2 ^ ceil(log2(spec_states.m/2));
  spec_states.dec           = 4;
  spec_states.sample_memory = zeros(1, spec_states.m);
  spec_states.lower_clip_dB = -30;

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
    Nmedian = 25;
    foff_median=zeros(1,Nmedian);

    % hf sim states
    
    path2 = zeros(1,hf_delay+M);
    sum_sig   = 0;
    sum_noise = 0;

    % state machine
    state = 0;
    fest_current = 0;
    fdelta = 5;
    candidate_thresh = 10;

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

      % optional HF channel sim

      if hf_sim
        path1 = tx_fdm .* conj(spread(f*M+1:f*M+M)');

        path2(1:hf_delay) = path2(M+1:hf_delay+M);
        path2(hf_delay+1:hf_delay+M) = tx_fdm .* conj(spread_2ms(f*M+1:f*M+M)');

        tx_fdm = hf_gain*(path1 + path2(1:M));
      end
      sum_sig += tx_fdm * tx_fdm';
      
      rx_fdm = real(tx_fdm);

      % AWGN noise

      noise = Nsd*randn(1,M); 
      sum_noise += noise * noise';
      rx_fdm += noise; 
      rx_fdm_log = [rx_fdm_log rx_fdm];

      % Fixed Delay

      rx_fdm_delay(1:Ndelay) = rx_fdm_delay(M+1:M+Ndelay);
      rx_fdm_delay(Ndelay+1:M+Ndelay) = rx_fdm; 

      % ------------------- Freq Offset Est -------------------

      % frequency offset estimation and correction, need to call
      % rx_est_freq_offset even in track mode to keep states updated

      [pilot prev_pilot pilot_lut_index prev_pilot_lut_index] = ...
          get_pilot(pilot_lut_index, prev_pilot_lut_index, M); 
      [foff_est S1 S2] = rx_est_freq_offset(rx_fdm_delay, pilot, prev_pilot, M);
      pilot_lpf1_log = [pilot_lpf1_log pilot_lpf1(Npilotlpf-M+1:Npilotlpf)];
      S1_log(f,:) = fftshift(S1);
     
      % raw estimate
      
      foff_log(ne,f) = foff_est;
      maxcf_log(ne,f) = max(abs(S1));

      % median filter post-processed

      foff_median(1:Nmedian-1) = foff_median(2:Nmedian);
      foff_median(Nmedian) = foff_est;
      foff_median_log(ne,f) = foff_coarse = median(foff_median);

      % state machine post-processed

      next_state = state;
      if state == 0
        if abs(foff_est - fest_current) > fdelta
          fest_candidate = foff_est;
          candidate_count = 0; 
          next_state = 1;
        end
      end
      if state == 1
         if abs(foff_est - fest_candidate) > fdelta
           next_state = 0;
         end
         candidate_count++; 
         if candidate_count > candidate_thresh
           fest_current = fest_candidate;
           next_state = 1;
        end
      end
      state = next_state;
      foff_statemach_log(ne,f) = fest_current;

      if (f > startup_delay) && (abs(foff_est - foff_hz) < allowable_error)
        hits++;
      end

      if length(EbNovec) == 1
        [spectrogram(f,:) spec_states] = spectrogram_update(rx_fdm, M, spec_states);
      end
    end

    % results for this EbNo value

    sim_out.foff_sd(ne) = std(foff_log(ne,startup_delay:frames));
    sim_out.hits = hits;
    sim_out.hits_percent = 100*sim_out.hits/(frames-startup_delay);
    sim_out.SNRvec(ne) = SNR;
    sim_out.tx_fdm_log = tx_fdm_log;
    sim_out.rx_fdm_log = rx_fdm_log;

    % noise we have measures is 4000 Hz wide, we want noise in 3000 Hz BW

    snr_meas = 10*log10(sum_sig/(sum_noise*4000/3000));

    printf("EbNo (dB): %5.2f  SNR: % -4.2f % -4.2f std dev (Hz): %3.2f  Hits: %d (%3.2f%%)\n", ...
           EbNo_dB, SNR, snr_meas, sim_out.foff_sd(ne), sim_out.hits, sim_out.hits_percent);

    % plots if single dimension vector

    if length(EbNovec) == 1
      figure(1)
      clf;
      subplot(311)
      plot(foff_log(ne,:))
      axis([1 frames -200 200]);
      ylabel("Freq offset estimate")
      subplot(312)
      plot(foff_median_log(ne,:))
      axis([1 frames -200 200]);
      ylabel("Freq offset estimate")
      subplot(313)
      plot(foff_statemach_log(ne,:))
      axis([1 frames -200 200]);
      xlabel("Frames")
      ylabel("Freq offset estimate")
      grid;

      figure(2)
      clf;
      plot(maxcf_log(ne,:))
      axis([1 frames -0 300]);
      xlabel("Frames")
      ylabel("max(abs(S1))")
      grid;

      figure(3)
      [n m] = size(S1_log);
      mesh(-200+400*(0:m-1)/256,1:n,abs(S1_log(:,:)));
      xlabel('Freq (Hz)'); ylabel('Frame num'); zlabel("max(abs(S1))")

      figure(4)
      clf
      [n m] = size(spectrogram);
      lower = floor(500*m/4000); upper = floor(2500*m/4000);
      mesh(lower:upper,1:n,spectrogram(:,lower:upper));
      xlabel('Freq (Hz)'); ylabel('Frame num'); zlabel('Amplitude (dB)');

      sim_out.spec = spectrogram;
      sim_out.tx_fdm_log = spectrogram;
    end
  end
end

% ---------------------------------------------------------------------
% Run Automated Tests
% ---------------------------------------------------------------------

figure(5); h=freqz(pilot_coeff,1,4000); plot(20*log10(abs(h(1:1000)))); grid

more off;

sim_in.EbNovec = 3;
sim_in.hf_sim = 1;
sim_in.hf_delay_ms = 2;
sim_in.delay = M/2;
sim_in.frames = Rs*10;
sim_in.foff_hz = 50;
sim_in.startup_delay = 10;
sim_in.allowable_error = 5;

sim_out = freq_off_est_test(sim_in);

% Test 1 - range of Eb/No (SNRs) in AWGN channel

function test1
  sim_in.EbNovec = 0:10;
  sim_in.delay = M/2;
  sim_in.frames = 20;
  sim_in.foff_hz = 50;
  sim_in.startup_delay = 10;
  sim_in.allowable_error = 5;

  sim_out = freq_off_est_test(sim_in);

  figure(4)
  clf
  subplot(211)
  plot(sim_in.EbNovec,sim_out.foff_sd)
  hold on;
  plot(sim_in.EbNovec,sim_out.foff_sd,'+')
  hold off;
  xlabel("Eb/No (dB)")
  ylabel("Std Dev")
  axis([(min(sim_in.EbNovec)-1) (max(sim_in.EbNovec)+1) -1 10]);

  subplot(212)
  plot(sim_out.SNRvec,sim_out.foff_sd)
  hold on;
  plot(sim_out.SNRvec,sim_out.foff_sd,'+')
  hold off;
  xlabel("SNR (dB)")
  ylabel("Std Dev")
  axis([(min(sim_out.SNRvec)-1)  (max(sim_out.SNRvec)+1) -1 10]);
end

% Test 2 - range of Eb/No (SNRs) in multipath channel

function test2
  sim_in.EbNovec = 0:10;
  sim_in.delay = 2;
  sim_in.hf_sim = 0;
  sim_in.hf_delay_ms = 2;
  sim_in.frames = Rs*10;
  sim_in.foff_hz = 0;
  sim_in.startup_delay = 10;
  sim_in.allowable_error = 5;

  sim_out = freq_off_est_test(sim_in);

  figure(5)
  clf
  subplot(211)
  plot(sim_in.EbNovec,sim_out.foff_sd)
  hold on;
  plot(sim_in.EbNovec,sim_out.foff_sd,'+')
  hold off;
  xlabel("Eb/No (dB)")
  ylabel("Std Dev")
  axis([(min(sim_in.EbNovec)-1) (max(sim_in.EbNovec)+1) -1 10]);
end

