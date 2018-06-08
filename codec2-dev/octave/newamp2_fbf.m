% newamp2_fbf.m
%
% Copyright David Rowe 2018
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Interactive Octave script to explore frame by frame operation of new amplitude
% modelling model version 2.
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README for
%    instructions.
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
%   $ cd ~/codec2-dev/octave
%   octave:14> newamp2_fbf("../build_linux/src/hts1a",50)


function newamp2_fbf(samname, f=73, varargin)
  more off;
  newamp2;
  Fs = 8000; 
  mode = "mel"; K = 20; correct_rate_K_en = 0;

  % command line arguments
  
  ind = arg_exists(varargin, "phase");
  phase_en = 0;
  if ind
    phase_en = 1;
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
    phase_source = 'Am';
  end
  
  % Keyboard loop --------------------------------------------------------------

  quant_en = 0;
  k = ' ';
  do 
    fg = 1;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    figure(fg++); clf; plot(s); axis([1 length(s) -20000 20000]);

    Wo = model(f,1); L = model(f,2); Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    K = 20;
    [rate_K_vec rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model(f,:), K, Fs, 'lanc');
    if correct_rate_K_en
      [tmp_ AmdB_] = resample_rate_L(model(f,:), rate_K_vec, rate_K_sample_freqs_kHz, Fs, "lancmel");
      [rate_K_vec_corrected orig_error error nasty_error_log nasty_error_m_log] = correct_rate_K_vec(rate_K_vec, rate_K_sample_freqs_kHz, AmdB, AmdB_, K, Wo, L, Fs);
      rate_K_vec = rate_K_vec_corrected;
    end

    % plots ----------------------------------
  
    figure(fg++); clf;
    l = sprintf(";rate %d AmdB;g+-", L);
    plot((1:L)*Wo*4000/pi, AmdB, l);
    axis([1 4000 -20 80]);
    hold on;

    rate_K_vec_ = rate_K_vec; 
    D = dct(rate_K_vec);
    if quant_en
       %rate_K_vec_ = huffman_quantise_rate_K(rate_K_vec);
       D_ = 16*round(D/16);
       rate_K_vec_ = idct(D_);
    else
       D_ = D;       
    end
    stem(rate_K_sample_freqs_kHz*1000, rate_K_vec_);
    
    % And .... back to rate L
    
    [model_ AmdB_] = resample_rate_L(model(f,:), rate_K_vec_, rate_K_sample_freqs_kHz, Fs, "lancmel");
    
    AmdB_ = AmdB_(1:L);
    sdL = std(abs(AmdB - AmdB_));

    plot((1:L)*Wo*4000/pi, AmdB_,";AmdB bar;r+-");
    l = sprintf(";error sd %3.2f dB;bk+-", sdL);
    plot((1:L)*Wo*4000/pi, (AmdB - AmdB_), l);
    hold off;

    if phase_en

      % est phase using HT

      fft_enc = 512;
      if strcmp(phase_source,'Am')
        Am_phase = model(f,3:(L+2));
        phase_est = determine_phase(model(f,:), 1, fft_enc);
      else
        Am_phase = model_(3:(L+2));
        phase_est = determine_phase(model_, 1, fft_enc);
      end
      Am_phase_dB = 20*log10(Am_phase);
      phase0 = zeros(1,L);
      for m=1:L
        b = round(m*Wo*fft_enc/(2*pi));
        phase0(m) = phase_est(b);
      end
      
      % plot amplitudes and phase for first 1kHz
      
      figure(fg++); clf;
      subplot(211);
      l = sprintf("+-;%s;",phase_source);
      plot((1:L)*Wo*4000/pi, Am_phase_dB(1:L),l);
      subplot(212);
      plot((1:L)*Wo*4000/pi, phase(f,1:L),'+-;orig;');
      hold on;
      plot((1:L)*Wo*4000/pi, phase0(1:L),'r+-;synth;');
      hold off;
      
      % simple synthesis using sinusoidal parameters

      figure(fg++); clf;
      N = 320;
      s = s_phase0 = zeros(1,N);
      for m=1:L
        s = s + Am_phase(m)*cos(m*Wo*(1:N) + phase(f,m));
        s_phase0 = s_phase0 + Am_phase(m)*cos(m*Wo*(1:N) + phase0(m));
      end
      subplot(211); plot(s); subplot(212); plot(s_phase0,'g');
    end

      figure(fg++);
      stem(D_)
      axis([0 K -50 50])
    
    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  u-qUant  c-Correct s-phSrc q-quit", f);
    fflush(stdout);
    k = kbhit();

    if k == 'n'
      f = f + 1;
    endif
    if k == 'b'
      f = f - 1;
    endif
    if k == 'u'
      if quant_en==0, quant_en=1, else quant_en=0;,end;
    end
    if k == 's'
      if strcmp(phase_source, 'Am')
        phase_source = 'Am bar'
      else
        phase_source = 'Am';
      end
    end
    if k == 'c'
      if correct_rate_K_en==0, correct_rate_K_en=1;, else correct_rate_K_en=0;,end;
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


