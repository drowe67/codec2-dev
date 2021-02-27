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

% Build a vector of Tx bursts in noise, one burst occurs every padded_burst_len samples

function [rx tx_preamble tx_postamble burst_len padded_burst_len ct_targets states] = generate_bursts(sim_in)
  config = ofdm_init_mode(sim_in.mode);
  states = ofdm_init(config);
  ofdm_load_const;

  tx_preamble = states.tx_preamble; tx_postamble = states.tx_postamble;
   
  Nbursts = sim_in.Nbursts;
  tx_bits = create_ldpc_test_frame(states, coded_frame=0);
  tx_burst = [tx_preamble ofdm_mod(states, tx_bits) tx_postamble];
  burst_len = length(tx_burst);
  tx_burst = ofdm_hilbert_clipper(states, tx_burst, tx_clip_en=0);
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
  mark_space_SNR_offset = 10*log10(burst_len/padded_burst_len);
  SNRdB_setpoint = sim_in.SNR3kdB + mark_space_SNR_offset;
  %printf("SNR3kdB: %f Burst offset: %f\n", sim_in.SNR3kdB, mark_space_SNR_offset)
  rx = channel_simulate(Fs, SNRdB_setpoint, sim_in.foff_Hz, sim_in.channel, tx);
endfunction


% Run an acquisition test, returning vectors of estimation errors.  Each tests runs on a segment of
% simulated signal that is known to contain a single burst.  We look for the peak of the acquisition
% metric in that segment.  This tests the core metric but is not (quite) a complete acquisition system as it
% doesn't operate frame by frame and operates on the assumption that the samples do contain a burst.
 
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
function ok = both_ok(dct, ttol_samples, dfoff, ftol_hz)
  Ntests = length(dct);
  ok = 0;
  for i = 1:Ntests
    if ((abs(dct(i)) < ttol_samples) && (abs(dfoff(i)) < ftol_hz))
      ok+=1;
    end
  end
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
  Pa = both_ok(dct, ttol_samples, dfoff, ftol_hz)/Ntests;
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


function results = acquisition_detector(states, rx, n, known_sequence)
  ofdm_load_const;
    
  % initial search over coarse grid
  tstep = 4; fstep = 5;
  [ct_est foff_est timing_mx] = est_timing_and_freq(states, rx(n:n+2*Nsamperframe-1), known_sequence, 
                                                    tstep, fmin = -50, fmax = 50, fstep);
  % refine estimate over finer grid                             
  fmin = foff_est - ceil(fstep/2); fmax = foff_est + ceil(fstep/2); 
  st = max(1, n + ct_est - tstep/2); en = st + Nsamperframe + tstep - 1;
  [ct_est foff_est timing_mx] = est_timing_and_freq(states, rx(st:en), known_sequence, 1, fmin, fmax, 1);
  ct_est += st;
  results.ct_est = ct_est; results.foff_est = foff_est; results.timing_mx = timing_mx;
end


function results = evaluate_candidate(states, det, i, Nsamperburstpadded, ct_target, foff_Hz, ttol_samples, ftol_hz)
  results.candidate = 0;
  if det.timing_mx > states.timing_mx_thresh
    % OK we have located a candidate peak
    
    % re-base ct_est to be wrt start of current burst reference frame
    ct_est = det.ct_est - (i-1)*Nsamperburstpadded;
        
    delta_ct = abs(ct_est-ct_target);
    delta_foff = det.foff_est-foff_Hz;
	
    ok = (abs(delta_ct) < ttol_samples) && (abs(delta_foff) < ftol_hz);
    
    results.candidate = 1; results.ct_est = ct_est; results.delta_ct = delta_ct; results.delta_foff = delta_foff; results.ok = ok;
  end    
endfunction


% test frame by frame acquisition algorithm

function Pa = frame_by_frame_acquisition_test(mode="datac1", Ntests=10, channel="awgn", SNR3kdB=100, foff_Hz=0, verbose_top=0) 
  sim_in.SNR3kdB = SNR3kdB;
  sim_in.channel = channel;
  sim_in.foff_Hz = foff_Hz;  
  sim_in.mode = mode;
  sim_in.Nbursts = Ntests;
  [rx tx_preamble tx_postamble Nsamperburst Nsamperburstpadded ct_targets states] = generate_bursts(sim_in);
  states.verbose = bitand(verbose_top,3);
  ofdm_load_const;
  
  timing_mx_log = []; ct_log = []; delta_ct_log = []; delta_foff_log = [];

  % allowable tolerance for acquistion
  ftol_hz = 2;              % we can sync up on this (todo: make mode selectable)
  ttol_samples = 0.006*Fs;  % CP length (todo: make mode selectable)
  target_acq = zeros(1,Ntests);
  
  for n=1:Nsamperframe:length(rx)-2*Nsamperframe
     pre = acquisition_detector(states, rx, n, tx_preamble);
     post = acquisition_detector(states, rx, n, tx_postamble);
     timing_mx_log = [timing_mx_log [pre.timing_mx; post.timing_mx]];
     
     % work out what burst we are evaluating
     i = ceil(n/Nsamperburstpadded); % i-th burst we are evaluating
     w = (i-1)*Nsamperburstpadded;   % offset of burst in s() for plotting purposes
     ct_target_pre = ct_targets(i);
     ct_target_post = ct_targets(i) + Nsamperburst - length(tx_preamble);

     pre_eval = evaluate_candidate(states, pre, i, Nsamperburstpadded, ct_target_pre, foff_Hz, ttol_samples, ftol_hz);
     post_eval = evaluate_candidate(states, post, i, Nsamperburstpadded, ct_target_post, foff_Hz, ttol_samples, ftol_hz);
     
     if pre_eval.candidate
       if pre_eval.ok == 0
          target_acq(i) = -1; % flag bad candidate
       end
       if pre_eval.ok && (target_acq(i) == 0)
          target_acq(i) = 1;  % flag a sucessful acquisition
       end
       delta_ct_log = [delta_ct_log pre_eval.delta_ct];
       delta_foff_log = [delta_foff_log pre_eval.delta_foff];
       ct_log = [ct_log w+pre_eval.ct_est];
       if states.verbose
         printf("Pre  i: %2d ct_est: %6d delta_ct: %6d foff_est: %5.1f timing_mx: %3.2f Acq: %d\n",
                 i, pre_eval.ct_est, pre_eval.delta_ct, pre.foff_est, pre.timing_mx, target_acq(i));
       end
     end  
     if post_eval.candidate
       if post_eval.ok == 0
          target_acq(i) = -1; % flag bad candidate
       end
       if post_eval.ok && (target_acq(i) == 0)
          target_acq(i) = 1;  % flag a sucessful acquisition
       end
       delta_ct_log = [delta_ct_log post_eval.delta_ct];
       delta_foff_log = [delta_foff_log post_eval.delta_foff];
       ct_log = [ct_log w+post_eval.ct_est];
       if states.verbose
         printf("Post i: %2d ct_est: %6d delta_ct: %6d foff_est: %5.1f timing_mx: %3.2f Acq: %d\n",
                i, post_eval.ct_est, post_eval.delta_ct, post.foff_est, post.timing_mx, target_acq(i));
       end
     end
            
  end
  
  if bitand(verbose_top,8)
    figure(1); clf; plot(timing_mx_log(1,:),'+-;preamble;'); 
               hold on; plot(timing_mx_log(2,:),'o-;postamble;'); hold off;
               title('mx log'); axis([0 length(timing_mx_log) 0 0.5]); grid;
    figure(4); clf; plot(real(rx)); axis([0 length(rx) -3E4 3E4]);
               hold on;
               plot(ct_log,zeros(1,length(ct_log)),'r+','markersize', 25, 'linewidth', 2);
               hold off; 
    figure(5); clf; plot_specgram(rx);
  end
  
  Pa = length(find(target_acq == 1))/Ntests;
  printf("%s %s SNR: %3.1f foff: %3.1f P(acq) = %3.2f\n", mode, channel, SNR3kdB, foff_Hz, Pa);
endfunction


% test frame by frame across modes, channels, and SNR (don't worry about sweeping freq)

function acquistion_curves_frame_by_frame_modes_channels_snr(Ntests=5, quick_test=0)
  modes={'datac0', 'datac1', 'datac3'};
  if quick_test
    Ntests = 5;
    channels={'awgn','mpp'}; SNR = [0 5];
  else
    channels={'awgn', 'mpm', 'mpp', 'notch'};
    SNR = [ -10 -5 -3.5 -1.5 0 1.5 3.5 5 7.5 10 15];
  end
  
  cc = ['b' 'g' 'k' 'c' 'm' 'r'];
  pt = ['+' '*' 'x' 'o' '+' '*'];
 
  for i=1:length(modes)
    figure(i); clf; hold on; title(sprintf("%s P(acquisition)", modes{i}));
  end
  
  for m=1:length(modes)
    figure(m);
    for c=1:length(channels)
      Pa_log = [];
      for s=1:length(SNR)
        Pa = frame_by_frame_acquisition_test(modes{m}, Ntests, channels{c}, SNR(s), foff_hz=0, verbose=1);
        Pa_log = [Pa_log Pa];
      end
      l = sprintf('%c%c-;%s;', cc(c), pt(c), channels{c}); 
      plot(SNR, Pa_log, l, 'markersize', 10);
    end
  end
  
  for i=1:length(modes)
    figure(i); grid;
    xlabel('SNR3k dB'); legend('location', 'southeast'); 
    xlim([min(SNR)-2 max(SNR)+2]); ylim([0 1.1]);
    print('-dpng', sprintf("%s_ofdm_dev_acq_curves_fbf_%s.png", datestr(clock(),"yyyy-mm-dd"), modes{i}));
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

%acquisition_test("datac3", Ntests=10, 'mpp', SNR3kdB=0, foff_hz=0, verbose=1+8);
%acquisition_histograms(mode="datac2", Ntests=3, channel='mpm', SNR3kdB=-5, foff=37, verbose=1+16)
%sync_metrics('freq')
%acquistion_curves("datac3", "mpp", Ntests=10)
%acquistion_curves_modes_channels(Ntests=25)
frame_by_frame_acquisition_test("datac1", Ntests=5, 'mpp', SNR3kdB=5, foff_hz=0, verbose=1+8);
%acquistion_curves_frame_by_frame_modes_channels_snr(Ntests=50, quick_test=0)
