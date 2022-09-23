% ratek2_fbf.m
%
% David Rowe 2022
%
% Rate K Experiment 2 - Filtering Am, resampling rate L<->K
%                     - interactive Octave script to explore frame by frame
%                       operation of rate K resampling
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README.md for
%    instructions.
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
%   $ cd ~/codec2-dev/octave
%   octave:14> ratek2_fbf("../build_linux/src/hts1a",50)


function ratek2_fbf(samname, f, resampler = 'spline')
  more off;

  newamp_700c;
  Fs = 8000; Nb = 20; K = 30; resampler = 'spline';

  % load up text files dumped from c2sim ---------------------------------------

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);
  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K);
  
  % Keyboard loop --------------------------------------------------------------

  k = ' '; energy = 1;
  do
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    figure(1); clf; plot(s); axis([1 length(s) -20000 20000]);

    Wo = model(f,1); F0 = Fs*Wo/(2*pi); L = model(f,2);
    Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    % Filter at rate L, Y = F(A)
    
    figure(2); clf;
    Y = zeros(1,L); YdB = zeros(1,L); hold on;
    for m=1:L
      h = generate_filter(m,F0,L,Nb);
      plot((1:L)*F0,h)
      if energy == 1
        Y(m) = sum(Am.^2 .* h);
        YdB(m) = 10*log10(Y(m));
      else
        YdB(m) = sum(AmdB .* h);
      end
    end
    hold off;

    % Resample to rate K, then back to rate L

    amodel = model(f,:); amodel(3:(L+2)) = 10 .^ (YdB/20);
    B = resample_const_rate_f(amodel, rate_K_sample_freqs_kHz, clip_en = 0, resampler);
    model_ = resample_rate_L(amodel, B, rate_K_sample_freqs_kHz, resampler);
    
    figure(3); clf;
    l = sprintf(";rate %d AmdB;g+-", L);
    plot((1:L)*Wo*4000/pi, AmdB, l);
    hold on;
    plot((1:L)*Wo*4000/pi, YdB, ';YdB;c+-');    
    stem(rate_K_sample_freqs_kHz*1000, B, ";rate K;b+-");
    Y_ = model_(1,3:(L+2)); YdB_ = 20*log10(Y_);
    Lmin = round(200/F0); Lmax = floor(3700/F0);
    E = sum((YdB(Lmin:Lmax) - YdB_(Lmin:Lmax)).^2)/(Lmax-Lmin+1);
    plot((1:L)*Wo*4000/pi, YdB_,";YdB hat;r+-");
    l = sprintf(";E %3.2f dB;bk+-", E);
    plot((Lmin:Lmax)*F0, (YdB(Lmin:Lmax) - YdB_(Lmin:Lmax)), l);
    axis([0 Fs/2 -10 80]);
    hold off;

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit p-png e-energy[%d]", f, energy);
    fflush(stdout);
    k = kbhit();

    if k == 'n'
      f = f + 1;
    endif
    if k == 'b'
      f = f - 1;
    endif
    if k == 'e'
      if energy == 1
        energy = 0;
      else
        energy = 1;
      end
    end
    if (k == 'p')
      [dir name ext]=fileparts("../build_linux/big_dog");
      set(gca, 'FontSize', 16);
      h = legend({"Rate L Am","Rate K Bm", "Rate L Am hat"}, "location", "north");
      legend("boxoff")
      set (h, "fontsize", 16);
      xlabel('Freq (Hz)'); ylabel('Amplitude (dB)');
      print(sprintf("ratek2_%s_%d",name,f),"-dpng","-S500,500");
    endif

  until (k == 'q')
  printf("\n");

endfunction
