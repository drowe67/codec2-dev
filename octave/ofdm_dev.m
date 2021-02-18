% ofdm_dev.m
% David Rowe Jan 2021
%
% Simulations used for development of HF data modem acquisition
%
% To run headless on a server:
%
%   DISPLAY=\"\" octave-cli --no-gui -qf ofdm_dev.m > 210218.txt &

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
function [rx tx_preamble burst_len padded_burst_len ct_targets states] = generate_bursts(sim_in)
  config = ofdm_init_mode(sim_in.mode);
  states = ofdm_init(config);
  ofdm_load_const;

  tx_preamble = generate_preamble(states);
  Nbursts = sim_in.Nbursts;
  tx_bits = create_ldpc_test_frame(states, coded_frame=0);
  tx_burst = [tx_preamble ofdm_mod(states, tx_bits) tx_preamble];
  burst_len = length(tx_burst);
  tx_burst = ofdm_hilbert_clipper(states, tx_burst, tx_clip_en=0);
  on_len = length(tx_burst);
  padded_burst_len = Fs+burst_len+Fs;
  
  tx = []; ct_targets = [];
  for f=1:Nbursts
    % 100ms of jitter in the burst start point
    jitter = floor(rand(1,1)*0.1*Fs);
    tx_burst_padded = [zeros(1,Fs+jitter) tx_burst zeros(1,Fs-jitter)];
    ct_targets = [ct_targets Fs+jitter];
    tx = [tx tx_burst_padded];
  end
  % adjust channel simulator SNR setpoint given (burst on length)/(sample length) ratio
  SNRdB_setpoint = sim_in.SNR3kdB + 10*log10(on_len/burst_len);
  rx = channel_simulate(Fs, SNRdB_setpoint, sim_in.foff_Hz, sim_in.channel, tx, verbose);
endfunction


% Run an acquisition test, returning vectors of estimation errors
function [delta_ct delta_foff timing_mx_log] = acquisition_test(mode="700D", Ntests=10, channel, SNR3kdB=100, foff_Hz=0, verbose_top=0)
  
  sim_in.SNR3kdB = SNR3kdB;
  sim_in.channel = channel;
  sim_in.foff_Hz = foff_Hz;  
  sim_in.mode = mode;
  sim_in.Nbursts = Ntests;
  [rx tx_preamble Nsamperburst Nsamperburstpadded ct_targets states] = generate_bursts(sim_in);
  states.verbose = bitand(verbose_top,3);
  ofdm_load_const;
  
  delta_ct = []; delta_foff = []; ct_log = []; timing_mx_log = []; 
   
  i = 1;
  states.foff_metric = 0;
  for w=1:Nsamperburstpadded:length(rx)
    [ct_est foff_est timing_mx] = est_timing_and_freq(states, rx(w:w+Nsamperburstpadded-1), tx_preamble, 
                                  tstep = 4, fmin = -50, fmax = 50, fstep = 5);
    fmin = foff_est-3; fmax = foff_est+3;
    st = w+ct_est; en = st + length(tx_preamble)-1; rx1 = rx(st:en);
    [tmp foff_est timing_mx] = est_timing_and_freq(states, rx1, tx_preamble, 
                                  tstep = 1, fmin, fmax, fstep = 1);

    % valid coarse timing could be pre-amble or post-amble
    ct_target1 = ct_targets(i);
    ct_target2 = ct_targets(i)+Nsamperburst-length(tx_preamble);
    %printf("  ct_target1: %d ct_target2: %d ct_est: %d\n", ct_target1, ct_target2, ct_est);
    ct_delta1 = ct_est-ct_target1;
    ct_delta2 = ct_est-ct_target2;
    adelta_ct = min([abs(ct_delta1) abs(ct_delta2)]);
    
    % log results
    delta_ct = [delta_ct adelta_ct];
    delta_foff = [delta_foff (foff_est-foff_Hz)];
    ct_log = [ct_log w+ct_est];
    timing_mx_log = [timing_mx_log; timing_mx];
    
    if states.verbose
      printf("i: %2d w: %8d ct_est: %6d delta_ct: %6d foff_est: %5.1f timing_mx: %3.2f\n",
              i++, w, ct_est, adelta_ct, foff_est, timing_mx);
    end

  end
  
  if bitand(verbose_top,8)
    figure(1); clf; plot(timing_mx_log,'+-'); title('mx log');
    figure(2); clf; plot(delta_ct,'+-'); title('delta ct');
    figure(3); clf; plot(delta_foff,'+-'); title('delta freq off');
    figure(4); clf; plot(real(rx)); hold on; plot(ct_log,zeros(1,length(ct_log)),'r+','markersize', 25, 'linewidth', 2); hold off;
    figure(5); clf; plot_specgram(rx);
  end
  
endfunction


% Helper function to count number of tests where both time and freq are OK
function P = both_ok(dct, ttol_samples, dfoff, ftol_hz)
  Ntests = length(dct);
  ok = 0;
  for i = 1:Ntests
    if ((abs(dct(i)) < ttol_samples) && (abs(dfoff(i)) < ftol_hz))
      ok+=1;
    end
  end
  P = ok/Ntests;
endfunction
  

#{
   Meausures aquisistion statistics for AWGN and HF channels
#}

function res = acquisition_histograms(mode="datac0", Ntests=10, channel = "awgn", SNR3kdB=100, foff=0, verbose=0)
  Fs = 8000;
  
  % allowable tolerance for acquistion

  ftol_hz = 2;              % we can sync up on this (todo: make mode selectable)
  ttol_samples = 0.006*Fs;  % CP length (todo: make mode selectable)

  [dct dfoff] = acquisition_test(mode, Ntests, channel, SNR3kdB, foff, verbose); 
  Pt = length(find (abs(dct) < ttol_samples))/length(dct);
  Pf = length(find (abs(dfoff) < ftol_hz))/length(dfoff);
  Pa = both_ok(dct, ttol_samples, dfoff, ftol_hz)
  printf("%s %s SNR: %3.1f foff: %3.1f P(time) = %3.2f  P(freq) = %3.2f P(acq) = %3.2f\n", mode, channel, SNR3kdB, foff, Pt, Pf, Pa);

  if bitand(verbose,16)
    figure(1); clf;
    hist(dct(find (abs(dct) < ttol_samples)))
    t = sprintf("Coarse Timing Error %s %s SNR = %3.2f foff = %3.1f", mode, channel, SNR3kdB, foff);
    title(t)
    figure(2); clf;
    hist(dfoff(find(abs(dfoff) < 2*ftol_hz)))
    t = sprintf("Coarse Freq Error %s %s SNR = %3.2f foff = %3.1f", mode, channel, SNR3kdB, foff);
    title(t);
  end

   res = [Pt Pf Pa];
endfunction


% plot some curves of Acquisition probability against EbNo and freq offset

function res_log = acquistion_curves(mode="datac1", channel='awgn', Ntests=10)
  SNR = [ -10 -5 0 5 10 15];
  foff = [-42 -7  0 49 ];
  %SNR = [-5 5]; foff = [ -2 2];
  cc = ['b' 'g' 'k' 'c' 'm' 'r'];
  pt = ['+' '*' 'x' 'o' '+' '*'];
  titles={'P(timing)', 'P(freq)', 'P(acq)'};
  png_suffixes={'time', 'freq', 'acq'};
  
  for i=1:length(titles)
    figure(i); clf; hold on; title(sprintf("%s %s %s",mode,channel,titles{i}));
  end

  res_log = []; % keep record of all results, one row per test, length(SNR) rows per freq step
  for f = 1:length(foff)
    afoff = foff(f);
    for e = 1:length(SNR)
      aSNR = SNR(e);
      res = acquisition_histograms(mode, Ntests, channel, aSNR, afoff, verbose=1);
      res_log = [res_log; res];
    end
    st = (f-1)*length(SNR)+1; en = st + length(SNR) - 1;
    for i=1:length(titles)
      figure(i); l = sprintf('%c%c-;%3.1f Hz;', cc(f), pt(f), afoff); plot(SNR, res_log(st:en,i), l, 'markersize', 10);
    end
  end
  
  for i=1:length(titles)
    figure(i); grid;
    xlabel('SNR3k dB'); legend('location', 'southeast'); 
    xlim([min(SNR)-2 max(SNR)+2]); ylim([0 1.1]);
    print('-dpng', sprintf("ofdm_dev_acq_curves_%s_%s_%s.png", mode, channel, png_suffixes{i}))
  end 
endfunction


% many acquisition curves across modes and channels

function res_log = acquistion_curves_modes_channels(Ntests=5)
  modes={'datac0', 'datac1', 'datac3'};
  channels={'awgn', 'mpm', 'mpp', 'notch'};
  res_log = [];
  for m=1:length(modes)
    for c=1:length(channels)
      res = acquistion_curves(modes{m}, channels{c}, Ntests);
      res_log = [res_log; res];
    end
  end
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


format;
more off;
pkg load signal;
graphics_toolkit ("gnuplot");
randn('seed',1);

% ---------------------------------------------------------
% choose simulation to run here 
% ---------------------------------------------------------

%acquisition_test("datac1", Ntests=5, 'notch', SNR3kdB=0, foff_hz=-38, verbose=1+8);
%acquisition_histograms(mode="datac2", Ntests=3, channel='mpm', SNR3kdB=-5, foff=37, verbose=1+16)
%sync_metrics('freq')
%acquistion_curves("datac0", "awgn", Ntests=3)
acquistion_curves_modes_channels(Ntests=25)
