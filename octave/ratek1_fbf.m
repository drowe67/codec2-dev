% ratek1_fbf.m
%
% David Rowe 2022
%
% Rate K Experiment 1 - L>K linear rateK resampling, interactive Octave script
% to explore frame by frame operation of rate K resampling
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README.md for
%    instructions.
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
%   $ cd ~/codec2-dev/octave
%   octave:14> ratek1_fbf("../build_linux/src/hts1a",50)


function ratek1_fbf(samname, f, resampler = 'spline')
  more off;

  newamp_700c;
  Fs = 8000;  K = 40;

  % load up text files dumped from c2sim ---------------------------------------

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);
  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);

  % pre-process
  [rate_K_surface sample_freqs_kHz] = resample_const_rate_f_lin(model(1:frames,:), K, resampler);

  % Keyboard loop --------------------------------------------------------------

  k = ' ';
  do
    fg = 1;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    figure(fg++); clf; plot(s); axis([1 length(s) -20000 20000]);

    Wo = model(f,1); F0 = Fs*Wo/(2*pi); L = model(f,2);
    Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    % plots ----------------------------------

    figure(fg++); clf;
    l = sprintf(";rate %d AmdB;g+-", L);
    plot((1:L)*Wo*4000/pi, AmdB, l);
    axis([1 4000 -20 80]);
    hold on;
    stem(sample_freqs_kHz*1000, rate_K_surface(f,:), ";rate K;b+-");

    % default
    rate_K_vec_ = rate_K_surface(f,:);

    % back to rate L
    model_(f,:) = resample_rate_L(model(f,:), rate_K_vec_, sample_freqs_kHz, resampler);
    Am_ = model_(f,3:(L+2)); AmdB_ = 20*log10(Am_);
    Lmin = round(200/F0); Lmax = floor(3700/F0);
    E = sum((AmdB(Lmin:Lmax) - AmdB_(Lmin:Lmax)).^2)/(Lmax-Lmin+1);

    plot((1:L)*Wo*4000/pi, AmdB_,";AmdB bar;r+-");
    l = sprintf(";E %3.2f dB;bk+-", E);
    plot((Lmin:Lmax)*F0, (AmdB(Lmin:Lmax) - AmdB_(Lmin:Lmax)), l);
    axis([0 Fs/2 -10 80]);
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
    if (k == 'p')
      [dir name ext]=fileparts("../build_linux/big_dog");
      set(gca, 'FontSize', 16);
      h = legend({"Rate L Am","Rate K Bm", "Rate L Am hat"}, "location", "north");
      legend("boxoff")
      set (h, "fontsize", 16);
      xlabel('Freq (Hz)'); ylabel('Amplitude (dB)');
      print(sprintf("ratek1_%s_%d",name,f),"-dpng","-S500,500");
    endif

  until (k == 'q')
  printf("\n");

endfunction
