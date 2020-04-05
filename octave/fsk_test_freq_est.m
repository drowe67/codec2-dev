% fsk_test_freq_est.m
% David Rowe April 2020

% test FSK frequency estimator, in particular at low Eb/No

fsk_lib;

function [f_log num_dud] = run_test(EbNodB = 10, num_frames=10)
  Fs = 8000;
  Rs = 100;
  M  = 4;
  bits_per_frame = 512;

  states = fsk_init(Fs,Rs,M);
  states.tx_real = 0;
  states.ftx = 900 + 2*states.Rs*(1:states.M);
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
  st = 1; f_log = [];
  for f=1:run_frames

    % extract nin samples from input stream
    nin = states.nin;
    en = st + states.nin - 1;

    % due to nin variations its possible to overrun buffer
    if en < length(rx)
      sf = rx(st:en);
      st += nin;
      states = est_freq(states, sf, states.M);
      f_log = [f_log; states.f];
    end
  end

  % Lets say that for a valid freq estimate, all four tones must be within 0.1*Rs of their tx freqeuncy
  num_dud = 0;
  for i=1:length(f_log)
    if sum(abs(f_log(i,:)-states.ftx) > 0.1*states.Rs)
      num_dud++;
    end
  end
end

function run_single
  EbNodB = 3;
  [f_log num_dud] = run_test(3, 100);

  percent_dud = 100*num_dud/length(f_log);
  printf("EbNodB: %4.2fdB tests: %d duds: %d %4.2f%% bad freq estimates\n", EbNodB, length(f_log), num_dud, percent_dud)

  figure(1); clf; plot(f_log)
  xlabel('Time (samples)'); ylabel('Frequency (Hz)');
  print("fsk_freq_est_single.png", "-dpng")
end

function run_curve

   EbNodB = 2:9;
   percent_log = [];
   for ne = 1:length(EbNodB)
      [f_log num_dud] = run_test(EbNodB(ne), 100);
      percent_dud = 100*num_dud/length(f_log);
      percent_log = [percent_log percent_dud];
      printf("EbNodB: %4.2f dB tests: %d duds: %3d %4.2f %% bad freq estimates\n", EbNodB(ne), length(f_log), num_dud, percent_dud)
  end
  
  figure(1); clf; plot(EbNodB, percent_log); grid;
  xlabel('Eb/No (dB)'); ylabel('% Errors');
  print("fsk_freq_est_curve.png", "-dpng")
end

graphics_toolkit("gnuplot");

% same results every time
%rand('state',1); 
%randn('state',1);

# choose one of these to run
#run_single
run_curve
