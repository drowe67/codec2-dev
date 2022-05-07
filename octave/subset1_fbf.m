
% subset1_fbf.m
%
% Interactive Octave script to explore frame by frame operation of newamp1
% subset spectral amplitude modelling.
%
% Usage:
%   $ c2sim big_dog_agc.s16 --dump big_dog --rateK \
%     --rateK_mean_min 0 --rateK_mean_max 60 --rateKout big_dog.f32
%   $ cd ~/codec2-dev/octave
%   octave:14> subset1_fbf("../script/big_dog",50)

function subset1_fbf(samname, f=73, varargin)
  more off;
  newamp_700c;
  graphics_toolkit("gnuplot");
  Fs = 8000;  K = 20;

  vq = 0; eq_en = 0; pf = 0;

  % load up text files dumped from c2sim ---------------------------------------

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);
  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);

  sample_freqs_kHz = mel_sample_freqs_kHz(K);
  rateKin_name = strcat(samname,".f32");
  rateKin = load_f32(rateKin_name,K);
  rateKout_name = strcat(samname,"_vq3.f32");
  rateKout = load_f32(rateKout_name,K);
  
  % Keyboard loop --------------------------------------------------------------

  k = ' ';
  do
    fg = 1;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    figure(fg++); clf; plot(s); axis([1 length(s) -20000 20000]);

    Wo = model(f,1); L = model(f,2); Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    % plots ----------------------------------

    figure(fg++); clf;
    l = sprintf(";rate %d AmdB;g+-", L);
    plot((1:L)*Wo*4000/pi, AmdB, l);
    axis([1 4000 -20 80]);
    hold on;

    rate_K_vec_ = rateKin(f,:);
    if vq
       rate_K_vec_ = rateKout(f,:);
    end
    stem(sample_freqs_kHz*1000, rate_K_vec_, ";rate K;b+-");
    
    % back to rate L
    model_(f,:) = resample_rate_L(model(f,:), rate_K_vec_, sample_freqs_kHz);
    Am_ = model_(f,3:(L+2)); AmdB_ = 20*log10(Am_);
    varL = var(AmdB - AmdB_);

    plot((1:L)*Wo*4000/pi, AmdB_,";AmdB bar;r+-");
    l = sprintf(";error var %3.2f dB;bk+-", varL);
    plot((1:L)*Wo*4000/pi, (AmdB - AmdB_), l);
    hold off;

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit  v-vq[%d] p-pf[%d] e-eq[%d]", f, vq, pf, eq_en);
    fflush(stdout);
    k = kbhit();

    if k == 'v'
       if vq == 0; vq = 1; else vq = 0; end
     endif
    if k == 'p'
       if pf == 0; pf = 1; else pf = 0; end
    endif
    if k == 'e'
       if eq_en == 0; eq_en = 1; else eq_en = 0; end
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
      if !ind && strcmp(v{i}, str)
        ind = i;
      end
    end
endfunction
