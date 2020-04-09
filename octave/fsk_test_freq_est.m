% fsk_test_freq_est.m
% David Rowe April 2020

% test FSK frequency estimator, in particular at low Eb/No

fsk_lib;

function [states f_log f_log2 num_dud1 num_dud2] = run_test(EbNodB = 10, num_frames=10)
  Fs = 8000;
  Rs = 100;
  M  = 4;
  bits_per_frame = 512;

  states = fsk_init(Fs,Rs,M);
  states.tx_real = 0;
  states.ftx = 900 + 2*states.Rs*(1:states.M);
  states.tx_tone_separation = 2*states.Rs;
  N = states.N;

  % Freq. estimator limits - keep these narrow to stop errors with low SNR 4FSK
  states.fest_fmin = 300;
  states.fest_fmax = 2200;
  states.fest_min_spacing = 100;

  EbNo = 10^(EbNodB/10);
  variance = states.Fs/(states.Rs*EbNo*states.bitspersymbol);

  tx_bits = round(rand(1,bits_per_frame*num_frames));
  tx = fsk_mod(states, tx_bits);
  noise = sqrt(variance/2)*randn(length(tx),1) + j*sqrt(variance/2)*randn(length(tx),1);
  rx = tx + noise;

  run_frames = floor(length(rx)/N)-1;
  st = 1; f_log = []; f_log2 = [];
  for f=1:run_frames

    % extract nin samples from input stream
    nin = states.nin;
    en = st + states.nin - 1;

    % due to nin variations it's possible to overrun buffer
    if en < length(rx)
      sf = rx(st:en);
      st += nin;
      states = est_freq(states, sf, states.M);
      f_log = [f_log; states.f]; f_log2 = [f_log2; states.f2];
    end
  end

  % Lets say that for a valid freq estimate, all four tones must be within 0.1*Rs of their tx freqeuncy
  num_dud1 = 0; num_dud2 = 0;
  for i=1:length(f_log)
    if sum(abs(f_log(i,:)-states.ftx) > 0.1*states.Rs)
      num_dud1++;
    end
    if sum(abs(f_log2(i,:)-states.ftx) > 0.1*states.Rs)
      num_dud2++;
    end
  end
end

function run_single(EbNodB = 3, num_frames = 10)
  [states f_log f_log2 num_dud1 num_dud2] = run_test(EbNodB, num_frames);
  
  percent_dud1 = 100*num_dud1/length(f_log);
  percent_dud2 = 100*num_dud2/length(f_log);
  printf("EbNodB: %4.2f dB tests: %3d duds1: %3d %5.2f %% duds1: %3d %5.2f %%\n",
         EbNodB, length(f_log), num_dud1, percent_dud1, num_dud2, percent_dud2)

  figure(1); clf;
  plot(f_log(:,1), 'linewidth', 2, 'b;peak;');
  hold on;
  plot(f_log(:,2:states.M), 'linewidth', 2, 'b');
  plot(f_log2(:,1),'linewidth', 2, 'r;mask;');
  plot(f_log2(:,2:states.M),'linewidth', 2, 'r');
  hold off;
  xlabel('Time (frames)'); ylabel('Frequency (Hz)');
  title(sprintf("EbNo = %4.2f dB", EbNodB));
  print("fsk_freq_est_single.png", "-dpng")
end

function run_curve

   EbNodB = 0:9;
   percent_log = [];
   for ne = 1:length(EbNodB)
      [states f_log f_log2 num_dud1 num_dud2] = run_test(EbNodB(ne), 100);
      percent_dud1 = 100*num_dud1/length(f_log);
      percent_dud2 = 100*num_dud2/length(f_log);
      percent_log = [percent_log; [percent_dud1 percent_dud2]];
      printf("EbNodB: %4.2f dB tests: %3d duds1: %3d %5.2f %% duds2: %3d %5.2f %% \n",
             EbNodB(ne), length(f_log), num_dud1, percent_dud1, num_dud2, percent_dud2)
  end
  
  figure(1); clf; plot(EbNodB, percent_log(:,1), 'linewidth', 2, '+-;peak;'); grid;
  hold on;  plot(EbNodB, percent_log(:,2), 'linewidth', 2, 'r+-;mask;'); hold off;
  xlabel('Eb/No (dB)'); ylabel('% Errors');
  print("fsk_freq_est_curve.png", "-dpng")
end

graphics_toolkit("gnuplot");
more off;

% same results every time
rand('state',1); 
randn('state',1);

# choose one of these to run
run_single(0)
run_curve
