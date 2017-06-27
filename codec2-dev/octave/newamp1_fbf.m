% newamp1_fbf.m
%
% Copyright David Rowe 2016
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Interactive Octave script to explore frame by frame operation of new amplitude
% modelling model.
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README for
%    instructions.
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
%   $ cd ~/codec2-dev/octave
%   octave:14> newamp1_fbf("../build_linux/src/hts1a",50)


function newamp1_fbf(samname, f=10, varargin)
  more off;

  newamp;
  melvq;

  Fs = 8000; rate_K_sample_freqs_kHz = [0.1:0.1:4]; K = length(rate_K_sample_freqs_kHz);
  quant_en = 0;

  % optional split VQ low freq quantiser

  ind = arg_exists(varargin, "vql");
  if ind
    quant_en = 1;
    vq_filename = varargin{ind+1};
    x = load(vq_filename); vq = x.vq;
    [vq_rows vq_cols] = size(vq); vq_st = 1; vq_en = vq_st + vq_cols - 1;
  end
  
  % optional split VQ high freq quantiser

  ind = arg_exists(varargin, "vqh");
  if ind
    quant_en = 1;
    vq_filename = varargin{ind+1};
    x = load(vq_filename); vq = x.vq;
    [vq_rows vq_cols] = size(vq); vq_st = 11; vq_en = K;
  end

  % load up text files dumped from c2sim ---------------------------------------

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);
  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);

  voicing_name = strcat(samname,"_pitche.txt");
  voicing = zeros(1,frames);
  
  if exist(voicing_name, "file") == 2
    pitche = load(voicing_name);
    voicing = pitche(:, 3);
  end

  % Keyboard loop --------------------------------------------------------------

  k = ' ';
  do 
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    figure(1); clf; plot(s); axis([1 length(s) -20000 20000]);

    Wo = model(f,1); L = model(f,2); Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    rate_K_vec = resample_const_rate_f(model(f,:), rate_K_sample_freqs_kHz, Fs);
    maskdB = determine_mask(rate_K_vec, rate_K_sample_freqs_kHz, rate_K_sample_freqs_kHz);
    rate_K_vec_no_mean = rate_K_vec - mean(rate_K_vec);

    % plots ----------------------------------
  
    figure(2); clf; 
    plot((1:L)*Wo*4000/pi, AmdB,";AmdB;g+-");
    axis([1 4000 -20 80]);
    hold on;
    plot(rate_K_sample_freqs_kHz*1000, rate_K_vec, ";rate K;b+-");
    %if quant_en != 2
    %  plot(rate_K_sample_freqs_kHz*1000, maskdB, ";mask;c+-");
    %end
    if quant_en == 2
      plot(rate_K_sample_freqs_kHz*1000, maskdB-rate_K_vec, ";maskdB - rate K;m+-");
    end
    #{
    masked_indx = find(maskdB >= rate_K_vec);
    unmasked_indx = find(maskdB < rate_K_vec);
    masked_indx = (vq_en-5:vq_en);
    unmasked_indx = (vq_st:vq_en-6);
    stem(rate_K_sample_freqs_kHz(unmasked_indx)*1000, 10*ones(1,length(unmasked_indx)), ";unmasked_index;ko-");
    #}

    if quant_en
      error = g = zeros(1, vq_rows);

      % find mse for each vector

      target = rate_K_vec_no_mean(vq_st:vq_en);
      weight_gain = 0.1;
      %masked_indx = find(maskdB(vq_en-5:vq_en) >= rate_K_vec(vq_en-5:vq_en));
      %masked_indx = find(masked_indx > (vq_cols-5))
      diff =  weights = diff_weighted = zeros(vq_rows, vq_cols);
      for i=1:vq_rows

        % work out gain for best match

        g(i) = sum(target - vq(i,:))/vq_cols;

        % Find weighted difference.  This allocated more importance
        % (error) to samples with higher energy, and stops really low
        % level harmonics from having any impact.  Note addition in dB
        % is multiplication in linear

        diff(i,:) = target - vq(i,:) - g(i);
        weights(i,:) = max(0.1,weight_gain .* (target + 20));
        
        diff_weighted(i,:) = diff(i,:) .* weights(i,:);

        error(i) = sum(abs(diff_weighted(i,:)));
      end

      [mn mn_ind] = min(error);

      rate_K_vec_no_mean_ = rate_K_vec_no_mean;
      rate_K_vec_no_mean_(vq_st:vq_en) = vq(mn_ind,:) + g(mn_ind);
      rate_K_vec_ = rate_K_vec_no_mean_ + mean(rate_K_vec);
      [model_ AmdB_] = resample_rate_L(model(f,:), rate_K_vec_, rate_K_sample_freqs_kHz, Fs);
      AmdB_ = AmdB_(1:L);
      if quant_en == 2
        plot(rate_K_sample_freqs_kHz(vq_st:vq_en)*1000, diff_weighted(mn_ind,:), ";diff;k+-");
      end
      plot((1:L)*Wo*4000/pi, AmdB_,";AmdB bar;r+-");
      hold off;

      figure (4); clf; plot(rate_K_sample_freqs_kHz(vq_st:vq_en)*1000, weights(mn_ind,:), ";weights;k+-");

      % sort and plot top m matches

      figure(3); clf;
      m = 4;
      [mse_list index_list] = sort(error);

      mse_list = mse_list(1:m);
      index_list = index_list(1:m);
      for i=1:m
        subplot(sqrt(m),sqrt(m),i);
        indx = index_list(i);
        plot(target + mean(rate_K_vec),'b+-');
        hold on;
        if index_list(i) == mn_ind
          plot(vq(indx,:) + g(indx) + mean(rate_K_vec),'r+-');
        else
          plot(vq(indx,:) + g(indx) + mean(rate_K_vec),'g+-');
        end
        plot(diff(indx,:), ";diff;k+-");
        plot(diff_weighted(indx,:), ";weighted diff;ko-");
        hold off;
        legend("location", "southwest");
      end
    end

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit  m-show_quant[%d]", f, quant_en);
    fflush(stdout);
    k = kbhit();

    if k == 'm'
      quant_en++;
      if quant_en > 2
        quant_en = 0;
      end
    endif
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
      if strcmp(v{i}, str)
        ind = i;
      end     
    end
endfunction
