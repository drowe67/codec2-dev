% twod_fbf.m
%
% Copyright David Rowe 2017
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Interactive Octave script to explore frame by frame operation with 2D error
% metrics
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README for
%    instructions.
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
%   $ cd ~/codec2-dev/octave
%   octave:14> twod_fbf("../build_linux/src/hts1a",50)


function twod_fbf(samname, f=73, varargin)
  more off;

  newamp;
  melvq;

  Fs = 8000; rate_K_sample_freqs_kHz = [0.1:0.1:4]; K = length(rate_K_sample_freqs_kHz);
  rate_K_sample_freqs_Hz = rate_K_sample_freqs_kHz*1000;

  % load up text files dumped from c2sim ---------------------------------------

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);
  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);
  
  % Keyboard loop --------------------------------------------------------------

  k = ' ';
  do 
    fg = 1;
    figure(fg++); clf;

    subplot(211,"position",[0.1 0.8 0.8 0.15])
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    plot(s);
    axis([1 length(s) -20000 20000]);

    Wo = model(f,1); L = model(f,2); Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi; Am_freqs_Hz = Am_freqs_kHz*1000;

    % Construct a "codebook" vector by inserting a single rate L
    % amplitude sample.  This will shift formants, showing up just the
    % the sort of situation we want our 2D error measure to work in.

    model_shifted = model(f,:);
    model_shifted(4:(L+2)) = model_shifted(3:(L+1));

    % resample to rate K

    rate_K_vec = resample_const_rate_f(model(f,:), rate_K_sample_freqs_kHz, Fs);
    rate_K_vec_ = resample_const_rate_f(model_shifted, rate_K_sample_freqs_kHz, Fs);

    % find closed point in rate K to rate L vector in terms of 2D distance

    twod_dist = twod_dist_f = twod_vec_x = twod_vec_y = zeros(1,L);
    weight_f = 3/100;
    for m=1:L
      min_dist = 1E32;
      for k=1:K
        dist  = (weight_f*(Am_freqs_Hz(m) - rate_K_sample_freqs_Hz(k))).^ 2;
        dist += (AmdB(m) - rate_K_vec_(k)).^2;
        if dist < min_dist
          min_dist = dist;
          min_k = k; 
        end
      end
      twod_dist_f(m) = Am_freqs_Hz(m);
      twod_dist(m) = min_dist;
      twod_vec_x(m) = rate_K_sample_freqs_Hz(min_k) - Am_freqs_Hz(m);
      twod_vec_y(m) = rate_K_vec_(min_k) - AmdB(m);
    end
    
    % plots ----------------------------------
  
    subplot(212,"position",[0.1 0.05 0.8 0.7])
    l = sprintf(";rate %d AmdB;g+-", L);
    plot(Am_freqs_Hz, AmdB, l);
    axis([1 4000 -20 80]);
    hold on;
    %plot(rate_K_sample_freqs_Hz, rate_K_vec, ";rate K;b+-");

    % And .... back to rate L
    
    [model_ AmdB_] = resample_rate_L(model(f,:), rate_K_vec_, rate_K_sample_freqs_kHz, Fs);
    AmdB_ = AmdB_(1:L);
    sdL = std(abs(AmdB - AmdB_));

    plot(Am_freqs_Hz, AmdB_,";AmdB bar;r+-");
    l = sprintf(";error sd %3.2f dB;bk+-", sdL);
    plot(Am_freqs_Hz, (AmdB - AmdB_), l);

    % 2D error and direction
    
    plot(Am_freqs_Hz, sqrt(twod_dist), ";2D error;c+-");
    for m=1:L
      plot([Am_freqs_Hz(m) Am_freqs_Hz(m)+twod_vec_x(m)], [AmdB(m) AmdB(m)+twod_vec_y(m)], 'c-');
      %[Am_freqs_Hz(m) Am_freqs_Hz(m)+twod_vec_x(m)]
    end
    hold off;

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit", f);
    fflush(stdout);
    k = kbhit();

    if k == 'n'
      f = f + 1;
    endif
    if k == 'b'
      f = f - 1;
    endif
  until (k == 'q')
  printf("\n");

endfunction

 
function ind = arg_exists(v, str) 
   ind = 0;
   for i=1:length(v)
      if !ind && strcmp(v{i}, str)
        ind = i;
      end     
    end
endfunction


