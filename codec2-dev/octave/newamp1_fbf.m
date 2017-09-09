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
  quant_en = 0; vq_search = "gain"; 
  mask_en = 0;
  nvq = 0; vq_start = [];
  quant_en = weight_en = 0;

  % specify one of more vqs, start index, vq file name, and search method

  ind = anind = arg_exists(varargin, "vq");
  while anind
    nvq++;
    vq_start = [vq_start varargin{ind+1}];
    avq_filename =  varargin{ind+2};
    if nvq < 2; vq_filename = avq_filename; else; vq_filename = [vq_filename; avq_filename]; end;
    avq_search =  varargin{ind+3};
    if nvq < 2; vq_search = avq_search; else; vq_search = [vq_search; avq_search]; end;
    printf("nvq %d vq_start: %d vq_filename: %s vq_search: %s\n", nvq, vq_start(nvq), avq_filename, avq_search);
    anind = arg_exists(varargin(ind+1:length(varargin)), "vq");
    if anind
      ind += anind;
    end
  end
  
  fit_order = 0;

  ind = arg_exists(varargin, "construct");
  if ind
    vq_search = varargin{ind};
  end
  ind = arg_exists(varargin, "construct_indep");
  if ind
    vq_search = varargin{ind};
  end
  
  % optional exploration of phase

  ind = arg_exists(varargin, "phase");
  phase_en = 0;
  if ind
    phase_en = 1;
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

  if phase_en
    phase_name = strcat(samname,"_phase.txt");
    phase = load(phase_name);
  end
  
  % Keyboard loop --------------------------------------------------------------

  k = ' ';
  do 
    fg = 1;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    figure(fg++); clf; plot(s); axis([1 length(s) -20000 20000]);

    Wo = model(f,1); L = model(f,2); Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    % remove constant gain term

    rate_K_vec = resample_const_rate_f(model(f,:), rate_K_sample_freqs_kHz, Fs);
    if fit_order == 0
      slope = 0; meanf = mean(rate_K_vec); 
      rate_K_vec_fit = rate_K_vec - meanf;
    end

    % plots ----------------------------------
  
    figure(fg++); clf;
    l = sprintf(";rate %d AmdB;g+-", L);
    plot((1:L)*Wo*4000/pi, AmdB, l);
    axis([1 4000 -20 80]);
    hold on;
    plot(rate_K_sample_freqs_kHz*1000, rate_K_vec, ";rate K;b+-");

    % default to the ideal
    
    rate_K_vec_ = rate_K_vec_fit; 
 
    if mask_en && nvq
      % experimental masking stuff that I can't seem to get to work
      maskdB = determine_mask(rate_K_vec, rate_K_sample_freqs_kHz, rate_K_sample_freqs_kHz, bark_model=1);
      plot(rate_K_sample_freqs_kHz*1000, maskdB, ";mask dB;c+-");
    end

    if strcmp(vq_search, "construct")
      [idx contrib errors b] = vq_construct_mg(rate_K_vec);
      rate_K_vec_ = contrib;
    end

    if strcmp(vq_search, "construct_indep")
      [idx contrib errors b] = vq_construct_indep_mg(rate_K_vec);
      rate_K_vec_ = contrib;
    end

    if nvq

      for i=1:nvq
        avq_filename = char(cellstr(vq_filename)(i));
        avq_search = char(cellstr(vq_search)(i));
        x = load(avq_filename); vq = x.vq; [vq_rows vq_cols] = size(vq);
        vq_st = vq_start(i); vq_en = vq_st + vq_cols - 1;
        printf("\nsplit VQ: %d vq_filename: %s vq_search: %s vq_st: %d vq_en: %d nVec: %d\n", i, avq_filename, avq_search, vq_st, vq_en, vq_rows);
        target = rate_K_vec_fit(vq_st:vq_en);

        if strcmp(avq_search, "mse")
          [idx contrib errors test_ g mg sl] = vq_search_mse(vq, target);
          rate_K_surface_fit_(f, vq_st:vq_en) = contrib;
        end

        if strcmp(avq_search, "gain")
          [idx contrib errors b] = vq_search_gain(vq, target, weights);
        end

        if strcmp(avq_search, "max")
          [idx contrib errors b] = vq_search_max(vq, target);
        end

        if strcmp(avq_search, "sg")
          [idx contrib errors b] = vq_search_sg(vq, target);
        end

        if strcmp(avq_search, "mg")
          [idx contrib errors b] = vq_search_mg(vq, target);
        end

        if strcmp(avq_search, "slope")
          [idx contrib errors b_log] = vq_search_slope(vq, target, "closed_quant_slope");
          rate_K_surface_fit_(f, vq_st:vq_en) = contrib;
          printf(" mg: %3.2f sl: %3.2f g: %3.2f \n", b_log(1), b_log(2), b_log(3));
          if quant_en
            % set slope to 0
            contrib1 = contrib;
            contrib = b_log(1)*vq(idx,:) + b_log(3);
            rate_K_vec_(vq_en+1:K) -= b_log(2)*vq_cols;
          end
        end

        rate_K_vec_(vq_st:vq_en) = contrib;
      end
      
      plot(rate_K_sample_freqs_kHz(vq_st:vq_en)*1000, contrib, 'm+-');
      if strcmp(vq_search, "para")
        plot(rate_K_sample_freqs_kHz(vq_st:vq_en)*1000, para_target, 'c+-');
        if quant_en
          plot(rate_K_sample_freqs_kHz(vq_st:vq_en)*1000, para, 'r+-');
        end
      end
      l = sprintf(";diff vq sd = %3.2f;k+-", std(target - contrib));
      plot(rate_K_sample_freqs_kHz(vq_st:vq_en)*1000, target - contrib, l);
    end

    % And .... back to rate L
    
    if (strcmp(vq_search, "construct") == 0) && (strcmp(vq_search, "construct_indep") == 0)
      rate_K_vec_ += meanf;
    end
    [model_ AmdB_] = resample_rate_L(model(f,:), rate_K_vec_, rate_K_sample_freqs_kHz, Fs);
    AmdB_ = AmdB_(1:L);
    sdL = std(abs(AmdB - AmdB_));

    plot((1:L)*Wo*4000/pi, AmdB_,";AmdB bar;r+-");
    if nvq == 0
      l = sprintf(";error sd %3.2f dB;bk+-", sdL);
      plot((1:L)*Wo*4000/pi, (AmdB - AmdB_), l);
    end
    hold off;

    if phase_en

      % est phase using HT

      Am_ = model(f,3:(L+2));
      fft_enc = 512;
      phase_est = determine_phase(model_, 1, fft_enc);
      phase0 = zeros(1,L);
      for m=1:L
        b = round(m*Wo*fft_enc/(2*pi));
        phase0(m) = phase_est(b);
      end
      
      % plot amplitudes and phase for first 1kHz
      
      figure(fg++); clf;
      subplot(211);
      plot((1:L)*Wo*4000/pi, AmdB_(1:L),'+-');
      subplot(212);
      plot((1:L)*Wo*4000/pi, phase(f,1:L),'+-');
      hold on;
      plot((1:L)*Wo*4000/pi, phase0(1:L),'r+-');
      hold off;
      
      % simple synthesis using sinusoidal parameters

      figure(fg++); clf;
      N = 320;
      s = s_phase0 = zeros(1,N);
      for m=1:L
        s = s + Am_(m)*cos(m*Wo*(1:N) + phase(f,m));
        s_phase0 = s_phase0 + Am_(m)*cos(m*Wo*(1:N) + phase0(m));
      end
      subplot(211); plot(s); subplot(212); plot(s_phase0,'g');
    end
    
    if quant_en
      figure(fg++); clf;
      plot(contrib1, 'b+-');
      hold on; plot(contrib,'r+'); hold off;
    end
    
    if weight_en
      figure(fg++); clf;
      subplot(211);
      plot((1:L)*Wo*4000/pi, AmdB,";AmdB;g+-");
      axis([1 4000 -20 80]);
      hold on;
      plot(rate_K_sample_freqs_kHz*1000, rate_K_vec, ";rate K;b+-");
      hold off;
      subplot(212);
      plot(rate_K_sample_freqs_kHz(vq_st:vq_en)*1000, weights);
      axis([1 4000 0 8]);
    end

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit  w-quant[%d]", f, quant_en);
    fflush(stdout);
    k = kbhit();

    if k == 'w'
      quant_en++;
      if quant_en == 2; quant_en = 0; end
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
      if !ind && strcmp(v{i}, str)
        ind = i;
      end     
    end
endfunction


