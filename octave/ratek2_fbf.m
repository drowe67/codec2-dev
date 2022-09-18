% ratek2_fbf.m
%
% David Rowe 2022
%
% Rate K Experiment 2 - Filtering Am, interactive Octave script
% to explore frame by frame operation of rate K resampling
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
  Fs = 8000; Nb = 20;

  % load up text files dumped from c2sim ---------------------------------------

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);
  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);

  % Keyboard loop --------------------------------------------------------------

  k = ' '; energy = 1;
  do
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    figure(1); clf; plot(s); axis([1 length(s) -20000 20000]);

    Wo = model(f,1); F0 = Fs*Wo/(2*pi); L = model(f,2);
    Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    % plots ----------------------------------

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
    
    figure(3); clf;
    l = sprintf(";rate %d AmdB;g+-", L);
    plot((1:L)*Wo*4000/pi, AmdB, l);
    axis([1 4000 -20 80]);
    hold on;
    plot((1:L)*Wo*4000/pi, YdB, ';YdB;r+-');    
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
