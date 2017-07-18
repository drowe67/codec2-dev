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


function newamp1_fbf(samname, f=73, varargin)
  more off;

  newamp;
  melvq;

  Fs = 8000; rate_K_sample_freqs_kHz = [0.1:0.1:4]; K = length(rate_K_sample_freqs_kHz);
  quant_en = 0; vq_search = "mse";

  % optional full band VQ

  ind = arg_exists(varargin, "vq");
  if ind
    quant_en = 1;
    vq_filename = varargin{ind+1};
    x = load(vq_filename); vq = x.vq;
    [vq_rows vq_cols] = size(vq); vq_st = 1; vq_en = vq_cols;
  end
  
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

  % different vq search algorithms

  ind = arg_exists(varargin, "vq_search");
  if ind
    vq_search = varargin{ind+1};
  end

  fit_order = 0;
  ind = arg_exists(varargin, "noslope");
  if ind
    fit_order = 1;
  end

  if quant_en
    printf("quant_en: %d vq_filename: %s vq_st: %d vq_en: %d vq_search: %s\n", 
           quant_en, vq_filename, vq_st, vq_en, vq_search);
  end

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
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    figure(1); clf; plot(s); axis([1 length(s) -20000 20000]);

    Wo = model(f,1); L = model(f,2); Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    % remove constant gain term

    rate_K_vec = resample_const_rate_f(model(f,:), rate_K_sample_freqs_kHz, Fs);
    if fit_order == 0
      slope = 0; b = mean(rate_K_vec); 
      rate_K_vec_fit = rate_K_vec - b;
    end

    % plots ----------------------------------
  
    figure(2); clf; 
    plot((1:L)*Wo*4000/pi, AmdB,";AmdB;g+-");
    axis([1 4000 -20 80]);
    hold on;
    plot(rate_K_sample_freqs_kHz*1000, rate_K_vec, ";rate K;b+-");

    if quant_en
      target = rate_K_vec_fit(vq_st:vq_en);
      
      if strcmp(vq_search, "mse")
        [idx contrib errors test_ g mg sl] = vq_search_mse(vq, target);
        rate_K_surface_fit_(f, vq_st:vq_en) = contrib;
      end

      if strcmp(vq_search, "slope")
        [idx contrib errors test_ g mg sl] = vq_search_slope(vq, target);
        rate_K_surface_fit_(f, vq_st:vq_en) = contrib;
      end

      rate_K_vec_ = rate_K_vec_fit; rate_K_vec_(vq_st:vq_en) = contrib;
      rate_K_vec_ += b;
      [model_ AmdB_] = resample_rate_L(model(f,:), rate_K_vec_, rate_K_sample_freqs_kHz, Fs);
      AmdB_ = AmdB_(1:L);

      sdL = std(AmdB - AmdB_);
      %printf("f: %d mn_ind: %d g: %3.2f sdK: %3.2f sdL: %3.2f\n", 
      %       f, mn_ind, g(mn_ind), error(mn_ind), sdL);

      plot(rate_K_sample_freqs_kHz(vq_st:vq_en)*1000, target - contrib, ";diff;k+-");
      plot((1:L)*Wo*4000/pi, AmdB_,";AmdB bar;r+-");
      hold off;

    end

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit  m-show_quant[%d]", f, quant_en);
    fflush(stdout);
    k = kbhit();

    if k == 'm'
      quant_en++;
      if quant_en == 2
        quant_en = 0;
      end
    endif
    if k == 'n'
      f = f + 1;
    endif
    if k == 'b'
      f = f - 1;
    endif
    if k == 'o'
      fit_order++;
      if fit_order == 2
        fit_order = 0;
      end
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


