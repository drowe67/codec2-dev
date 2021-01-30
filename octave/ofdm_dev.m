% ofdm_dev.m
% David Rowe Jan 2021
%
% Simulations used for development of HF data modem acquisition

ofdm_lib;
channel_lib;

% build a single modem frame preamble vector
function tx_preamble = generate_preamble(states)
  % tweak local copy of states so we can generate a 1 modem-frame packet
  states.Np = 1; states.Nbitsperpacket = states.Nbitsperframe;
  preamble_bits = rand(1,states.Nbitsperframe) > 0.5;
  tx_preamble = ofdm_mod(states, preamble_bits);
endfunction

% build a vector of Tx bursts in noise
function [rx tx_preamble states] = generate_bursts(sim_in)
  config = ofdm_init_mode(sim_in.mode);
  states = ofdm_init(config);
  ofdm_load_const;

  tx_preamble = generate_preamble(states);
  Nbursts = sim_in.Nbursts;
  tx_bits = create_ldpc_test_frame(states, coded_frame=0);
  tx = [tx_preamble ofdm_mod(states, tx_bits)];

  rx = [];
  for f=1:Nbursts
    arx = ofdm_clip_channel(states, tx, sim_in.SNR3kdB, sim_in.channel, sim_in.foff_Hz, 0);
    rx = [rx arx];
  end
endfunction

% Run an acquisition test, returning vectors of estimation errors
function [delta_ct delta_foff timing_mx_log] = acquisition_test(mode="700D", Ntests=10, SNR3kdB=100, foff_Hz=0, channel, verbose_top=0)
  
  sim_in.SNR3kdB = SNR3kdB;
  sim_in.channel = channel;
  sim_in.foff_Hz = foff_Hz;  
  sim_in.mode = mode;
  sim_in.Nbursts = Ntests;
  [rx tx_preamble states] = generate_bursts(sim_in);
  states.verbose = bitand(verbose_top,3);
  ofdm_load_const;
  
  delta_ct = []; delta_foff = []; ct_log = []; timing_mx_log = []; 
  
  Nsamperburst = Fs + Nsamperframe*(Np+1) + Fs/2; 
  printf("Nsamperburst: %d length: %d %d\n", Nsamperburst, length(rx), length(tx_preamble));
 
  i = 1;
  states.foff_metric = 0;
  for w=1:Nsamperburst:length(rx)
    [ct_est timing_valid timing_mx] = est_timing(states, rx(w:w+Nsamperburst-1), tx_preamble, 1, dual=0);
    foff_est = est_freq_offset_known_corr(states, rx(w:w+Nsamperburst-1), tx_preamble, ct_est, dual=0);
    if states.verbose
      printf("i: %2d w: %8d ct_est: %4d foff_est: %5.1f timing_mx: %3.2f timing_vld: %d\n", i++, w, ct_est, foff_est, timing_mx, timing_valid);
    end

    % valid coarse timing ests are modulo Nsamperframe

    ct_target = Fs+1;
    delta_ct = [delta_ct ct_est-ct_target];
    delta_foff = [delta_foff (foff_est-foff_Hz)];
    ct_log = [ct_log w+ct_est];
    timing_mx_log = [timing_mx_log; timing_mx];
  end
  
  if bitand(verbose_top,8)
    figure(1); clf; plot(timing_mx_log,'+-'); title('mx log');
    figure(2); clf; plot(delta_ct,'+-'); title('delta ct');
    figure(3); clf; plot(delta_foff,'+-'); title('delta freq off');
    figure(5); clf; plot(rx); hold on; plot(ct_log,zeros(1,length(ct_log)),'r+','markersize', 25, 'linewidth', 2); hold off;
  end
  
endfunction


#{

   Generates aquisistion statistics for AWGN and HF channels for
   continuous signals. Probability of acquistion is what matters,
   e.g. if it's 50% we can expect sync within 2 frames.

#}

function res = acquisition_histograms(mode="700D", fine_en = 0, foff, EbNoAWGN=-1, EbNoHF=3, verbose=1)
  Fs = 8000;
  Ntests = 100;

  % allowable tolerance for acquistion

  ftol_hz = 1.5;            % we can sync up on this
  ttol_samples = 0.002*Fs;  % 2ms, ie CP length

  % AWGN channel at uncoded Eb/No operating point

  [dct dfoff] = acquisition_test(mode, Ntests, EbNoAWGN, foff, 0, fine_en);

  % Probability of acquistion is what matters, e.g. if it's 50% we can
  % expect sync within 2 frames

  PtAWGN = length(find (abs(dct) < ttol_samples))/length(dct);
  printf("AWGN P(time offset acq) = %3.2f\n", PtAWGN);
  if fine_en == 0
    PfAWGN = length(find (abs(dfoff) < ftol_hz))/length(dfoff);
    printf("AWGN P(freq offset acq) = %3.2f\n", PfAWGN);
  end

  if verbose
    figure(1); clf;
    hist(dct(find (abs(dct) < ttol_samples)))
    t = sprintf("Coarse Timing Error AWGN EbNo = %3.2f foff = %3.1f", EbNoAWGN, foff);
    title(t)
    if fine_en == 0
      figure(2)
      hist(dfoff(find(abs(dfoff) < 2*ftol_hz)))
      t = sprintf("Coarse Freq Error AWGN EbNo = %3.2f foff = %3.1f", EbNoAWGN, foff);
      title(t);
    end
 end

  % HF channel at uncoded operating point

  [dct dfoff] = acquisition_test(mode, Ntests, EbNoHF, foff, 1, fine_en);

  PtHF = length(find (abs(dct) < ttol_samples))/length(dct);
  printf("HF P(time offset acq) = %3.2f\n", PtHF);
  if fine_en == 0
    PfHF = length(find (abs(dfoff) < ftol_hz))/length(dfoff)
    printf("HF P(freq offset acq) = %3.2f\n", PfHF);
  end

  if verbose
    figure(3); clf;
    hist(dct(find (abs(dct) < ttol_samples)))
    t = sprintf("Coarse Timing Error HF EbNo = %3.2f foff = %3.1f", EbNoHF, foff);
    title(t)
    if fine_en == 0
      figure(4)
      hist(dfoff(find(abs(dfoff) < 2*ftol_hz)))
      t = sprintf("Coarse Freq Error HF EbNo = %3.2f foff = %3.1f", EbNoHF, foff);
      title(t);
    end
  end
  
  res = [PtAWGN PfAWGN PtHF PfHF];
endfunction


% plot some curves of Acquisition probability against EbNo and freq offset

function acquistion_curves(mode="700D")

  EbNo = [-1 2 5 8 20];
  %foff = [-20 -15 -10 -5 0 5 10 15 20];
  foff = [-15 -5 0 5 15];
  cc = ['b' 'g' 'k' 'c' 'm'];
  
  figure(1); clf; hold on; title('P(timing) AWGN'); xlabel('Eb/No dB'); legend('location', 'southeast');
  figure(2); clf; hold on; title('P(freq) AWGN'); xlabel('Eb/No dB'); legend('location', 'southeast');
  figure(3); clf; hold on; title('P(timing) HF'); xlabel('Eb/No dB'); legend('location', 'southeast');
  figure(4); clf; hold on; title('P(freq) HF'); xlabel('Eb/No dB'); legend('location', 'southeast');

  for f=1:4
    ylim([0 1]);
  end
  
  for f = 1:length(foff)
    afoff = foff(f);
    res_log = [];
    for e = 1:length(EbNo)
      aEbNo = EbNo(e);
      res = zeros(1,4);
      res = acquisition_histograms(mode, fine_en = 0, afoff, aEbNo, aEbNo+4, verbose = 0);
      res_log = [res_log; res];
    end
    figure(1); l = sprintf('%c+-;%3.1f Hz;', cc(f), afoff); plot(EbNo, res_log(:,1), l);
    figure(2); l = sprintf('%c+-;%3.1f Hz;', cc(f), afoff); plot(EbNo, res_log(:,3), l);
    figure(3); l = sprintf('%c+-;%3.1f Hz;', cc(f), afoff); plot(EbNo+4, res_log(:,2), l);
    figure(4); l = sprintf('%c+-;%3.1f Hz;', cc(f), afoff); plot(EbNo+4, res_log(:,4), l);
  end
  
  figure(1); print('-deps', '-color', sprintf("ofdm_dev_acq_curves_time_awgn_%s.eps", mode))
  figure(2); print('-deps', '-color', sprintf("ofdm_dev_acq_curves_freq_awgn_%s.eps", mode))
  figure(3); print('-deps', '-color', sprintf("ofdm_dev_acq_curves_time_hf_%s.eps", mode))
  figure(4); print('-deps', '-color', sprintf("ofdm_dev_acq_curves_freq_hf_%s.eps", mode))
endfunction


% Used to develop sync state machine - in particular a metric to show
% we are out of sync, or have sync with a bad freq offset est, or have
% lost modem signal

function sync_metrics(mode = "700D", x_axis = 'EbNo')
  Fs      = 8000;
  Ntests  = 4;
  f_offHz = [-25:25];
  EbNodB  = [-10 0 3 6 10 20];
  %f_offHz = [-5:5:5];
  %EbNodB  = [-10 0 10];
  cc = ['b' 'g' 'k' 'c' 'm' 'b'];
  pt = ['+' '+' '+' '+' '+' 'o'];
    
  mean_mx1_log = mean_dfoff_log = [];
  for f = 1:length(f_offHz)
    af_offHz = f_offHz(f);
    mean_mx1_row = mean_dfoff_row = [];
    for e = 1:length(EbNodB)
      aEbNodB = EbNodB(e);
      [dct dfoff timing_mx_log] = acquisition_test(mode, Ntests, aEbNodB, af_offHz);
      mean_mx1 = mean(timing_mx_log(:,1));
      printf("f_offHz: %5.2f EbNodB: % 6.2f mx1: %3.2f\n", af_offHz, aEbNodB, mean_mx1);
      mean_mx1_row = [mean_mx1_row mean_mx1];
      mean_dfoff_row = [mean_dfoff_row mean(dfoff)];
    end
    mean_mx1_log = [mean_mx1_log; mean_mx1_row];
    mean_dfoff_log = [mean_dfoff_log; mean_dfoff_row];
  end

  figure(1); clf; hold on; grid;
  if strcmp(x_axis,'EbNo')
    for f = 1:length(f_offHz)
      if f == 2, hold on, end;
      leg1 = sprintf("b+-;mx1 %4.1f Hz;", f_offHz(f));
      plot(EbNodB, mean_mx1_log(f,:), leg1)
    end
    hold off;
    xlabel('Eb/No (dB)');
    ylabel('Coefficient')
    title('Pilot Correlation Metric against Eb/No for different Freq Offsets');
    legend("location", "northwest"); legend("boxoff");
    axis([min(EbNodB) max(EbNodB) 0 1.2])
    print('-deps', '-color', "ofdm_dev_pilot_correlation_ebno.eps")
  end

  if strcmp(x_axis,'freq')
    % x axis is freq

    for e = length(EbNodB):-1:1
      leg1 = sprintf("%c%c-;mx1 %3.0f dB;", cc(e), pt(e), EbNodB(e));
      plot(f_offHz, mean_mx1_log(:,e), leg1)
    end
    hold off;
    xlabel('freq offset (Hz)');
    ylabel('Coefficient')
    title('Pilot Correlation Metric against Freq Offset for different Eb/No dB');
    legend("location", "northwest"); legend("boxoff");
    axis([min(f_offHz) max(f_offHz) 0 1])
    print('-deps', '-color', "ofdm_dev_pilot_correlation_freq.eps")

    mean_dfoff_log
 
    figure(2); clf;
    for e = 1:length(EbNodB)
      if e == 2, hold on, end;
      leg1 = sprintf("+-;mx1 %3.0f dB;", EbNodB(e));
      plot(f_offHz, mean_dfoff_log(:,e), leg1)
    end
    hold off;
    xlabel('freq offset (Hz)');
    ylabel('Mean Freq Est Error')
    title('Freq Est Error against Freq Offset for different Eb/No dB');
    axis([min(f_offHz) max(f_offHz) -5 5])
  end
  
endfunction


% ---------------------------------------------------------
% choose simulation to run here 
% ---------------------------------------------------------

format;
more off;
pkg load signal;
randn('seed',1);

acquisition_test("datac1", Ntests=10, SNR3kdB=10, foff_hz=0, 'mpp', verbose=1+8);
%acquisition_histograms("700D", fin_en=0, foff_hz=-15, EbNoAWGN=-1, EbNoHF=3)
%sync_metrics('freq')
%acquisition_dev(Ntests=10, EbNodB=100, foff_hz=0)
%acquistion_curves
