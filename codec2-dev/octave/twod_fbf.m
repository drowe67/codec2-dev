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
  
  % Variety of tests to exercise new 2D error measure
  % default test1

  test_type = 1;
  ind = arg_exists(varargin, "test2");
  if ind
    test_type = 2;
    
    % take hts1a frame 73 and construct two VQ entries
    % i) HPF version ii) missing format, which means fix to one frame
    % show 1D and 2D error per vec

    % i) HPF: attenutate first few samples, attn as function of frequency, so 0-500Hz,
    % straight line, then flat after cut off.
    % y = mx + c, y = 0 at Xc, y = m(Xc)+c, m = (y - c)Xc = -c/Xc;
    
    c = -20; Xc= 500; m = -c/Xc;
    f = 73; Wo = model(f,1); L = model(f,2); Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    x = (1:L)*Fs*Wo/(2*pi);
    y = m*x + c;
    y(find(y>=0)) = 0;
    AmdB += y;
    
    % now resample to rate K and add to VQ table

    amodel = model(f,:); amodel(3:(L+2)) = 10 .^ (AmdB/20);
    rate_K_vec = resample_const_rate_f(amodel, rate_K_sample_freqs_kHz, Fs);
    vq = rate_K_vec;

    % ii) zero out a formant between 2000 and 2500 Hz

    Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    m_st = round(1900*2*pi/(Fs*Wo)); m_en = round(2200*2*pi/(Fs*Wo));
    AmdB(m_st:m_en) = AmdB(m_st-1);
    amodel = model(f,:); amodel(3:(L+2)) = 10 .^ (AmdB/20);
    rate_K_vec = resample_const_rate_f(amodel, rate_K_sample_freqs_kHz, Fs);
    vq = [vq; rate_K_vec];

    % Now search the two element vq, using the orginal as the target.
    % We would like vq(1,:) to be chosen, as a little HPF doesn't
    % affect the speech much, but losing a fomant does.  However we
    % expect the 1D model to choose vq(2,:), as the HPF distortion
    % will affect the metric moe than the missing formant.

    % need the 1D and 2D cost funcs in functions  Print chosen vectors, plot
    % 2D vectors
    
  end

  
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

    if test_type == 1
      % Construct a "codebook" vector by inserting a single rate L
      % amplitude sample.  This will shift formants, showing up just the
      % the sort of situation we want our 2D error measure to work in.

      model_shifted = model(f,:);
      model_shifted(4:(L+2)) = model_shifted(3:(L+1));

      % resample to rate K

      rate_K_vec = resample_const_rate_f(model(f,:), rate_K_sample_freqs_kHz, Fs);
      rate_K_vec_ = resample_const_rate_f(model_shifted, rate_K_sample_freqs_kHz, Fs);

      % find closed point in rate K to rate L vector in terms of 2D distance

      [twod_dist twod_dist_f twod_vec_x twod_vec_y ] = determine_twod_dist(AmdB, Am_freqs_Hz, rate_K_vec_, rate_K_sample_freqs_Hz);
    end
    
    % plots ----------------------------------
  
    subplot(212,"position",[0.1 0.05 0.8 0.7])
    l = sprintf(";rate L=%d AmdB;g+-", L);
    plot(Am_freqs_Hz, AmdB, l);
    axis([1 4000 -20 80]);
    hold on;

    if test_type == 1
      plot(rate_K_sample_freqs_Hz, rate_K_vec_, ";rate K v;b+-");

     % plot 2D error and direction

      plot(Am_freqs_Hz, sqrt(twod_dist), ";2D error;c+-");
      for m=1:L
        plot([Am_freqs_Hz(m) Am_freqs_Hz(m) - twod_vec_x(m)], [AmdB(m) AmdB(m) - twod_vec_y(m)], 'm-', 'linewidth', 2);
      end
    end
    if test_type == 2
      plot(rate_K_sample_freqs_Hz, vq(1,:), ";vq1;r+-");
      plot(rate_K_sample_freqs_Hz, vq(2,:), ";vq2;bk+-");

      rate_K_vec = resample_const_rate_f(model(f,:), rate_K_sample_freqs_kHz, Fs);
      
      % evaulate 2D cost function of target against two vectors
      
      [twod_dist1 twod_dist_f1 twod_vec_x twod_vec_y ] = determine_twod_dist(AmdB, Am_freqs_Hz, vq(1,:), rate_K_sample_freqs_Hz);
      [twod_dist2 twod_dist_f2 twod_vec_x twod_vec_y ] = determine_twod_dist(AmdB, Am_freqs_Hz, vq(2,:), rate_K_sample_freqs_Hz);

      figure(2);
      subplot(211)
      plot(rate_K_sample_freqs_Hz, abs(rate_K_vec - vq(1,:)),';vq1;b');
      hold on; plot(rate_K_sample_freqs_Hz, abs(rate_K_vec - vq(2,:)),';vq2;g'); hold off;
      title('1D Distance');
      
      subplot(212)
      plot(twod_dist_f1, twod_dist1,';vq1;b');
      hold on; plot(twod_dist_f2, twod_dist2,';vq2;g'); hold off;
      title('2D Distance');

      twoD1 = sum(twod_dist1);  twoD2 = sum(twod_dist2);
      oneD1 = sum(abs(rate_K_vec - vq(1,:))); oneD2 = sum(abs(rate_K_vec - vq(2,:)));
      [tmp, oneD_choice] = min([oneD1 oneD2]);
      [tmp, twoD_choice] = min([twoD1 twoD2]);
      
      printf("Vector   1D Distance   2D Distance\n");
      printf("----------------------------------\n");
      printf("vq(1,:):        %3.1f          %3.1f\n", oneD1, twoD1);
      printf("vq(2,:):        %3.1f          %3.1f\n", oneD2, twoD2);
      printf("Choice.:        %d                %d\n", oneD_choice, twoD_choice);
    end
    
    hold off;

    % interactive menu ------------------------------------------

    if test_type == 1
      printf("\rframe: %d  menu: n-next  b-back  q-quit", f);
      fflush(stdout);
      k = kbhit();

      if k == 'n'; f = f + 1; endif
      if k == 'b'; f = f - 1; endif
    else
      k = 'q';
    end
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


function [twod_dist twod_dist_f twod_vec_x twod_vec_y ] = determine_twod_dist(AmdB, Am_freqs_Hz, rate_K_vec_, rate_K_sample_freqs_Hz)

    L = length(AmdB); K = length(rate_K_vec_);
    
    twod_dist = twod_dist_f = twod_vec_x = twod_vec_y = zeros(1,L);
    weight_f = 0.05;

    for m=1:L

      % OK lets find two closest points to m-th point in rate L target

      dist  = (weight_f*(Am_freqs_Hz(m) - rate_K_sample_freqs_Hz)).^ 2;
      dist += (AmdB(m) - rate_K_vec_).^2;
      [tmp ind] = sort(dist);

      ind1 = ind(1);
      if ind(1) == 1
        ind2 = 2;
      elseif ind(1) == K
        ind1 = K-1; ind2 = K;
      else
        if dist(ind1-1) < dist(ind1+1)
          ind2 = ind1; ind1 = ind1-1;
        else
          ind2 = ind1+1;
        end
      end

      % construct x and p vectors

      %printf("m: %d ind11: %d ind2: %d\n", m, ind1, ind2);
      x_freq = weight_f*(rate_K_sample_freqs_Hz(ind2) - rate_K_sample_freqs_Hz(ind1));
      x_amp = rate_K_vec_(ind2) - rate_K_vec_(ind1);
      x = [x_freq x_amp];
      p_freq = weight_f*(Am_freqs_Hz(m) - rate_K_sample_freqs_Hz(ind1));
      p_amp = AmdB(m) - rate_K_vec_(ind1);
      p = [p_freq p_amp];

      % find gain g to make x orthogonal to p

      g = x*p'/(x*x'); e = p - g*x;
      
      twod_dist_f(m) = Am_freqs_Hz(m);
      twod_dist(m) = norm(e);
      twod_vec_x(m) = e(1)/weight_f;
      twod_vec_y(m) = e(2);
    end
endfunction
