% ratek2_nb_fbf.m
%
% David Rowe 2022
%
% Rate K Experiment 2 - Filtering Am, resampling rate L<->K
%                     - interactive Octave script to explore frame by frame
%                       operation of rate K resampling
%                     - this version used to plot Ym for different Nb
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README.md for
%    instructions.
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
%   $ cd ~/codec2-dev/octave
%   octave:14> ratek2_fbf("../build_linux/src/hts1a",50)


function ratek2_nb_fbf(samname, f, resampler = 'spline')
  more off;

  newamp_700c;
  Fs = 8000; Nb= [10 20];

  % load up text files dumped from c2sim ---------------------------------------

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);
  
  % Keyboard loop --------------------------------------------------------------

  k = ' ';
  do
    Wo = model(f,1); F0 = Fs*Wo/(2*pi); L = model(f,2);
    Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    % Filter at rate L, Y = F(A)
    
    Y = zeros(length(Nb),L); YdB = zeros(length(Nb),L);
    for m=1:L
      for i=1:length(Nb)
        h = generate_filter(m,F0,L,Nb(i));
        Y(i,m) = sum(Am.^2 .* h);
        YdB(i,m) = 10*log10(Y(i,m));
      end
    end
 
    figure(1); clf; l=1;
    leg{l} = sprintf("rate %d AmdB", L); l++;
    plot((1:L)*Wo*4000/pi, AmdB, '+-');
    hold on;
    for i=1:length(Nb)
      plot((1:L)*Wo*4000/pi, YdB(i,:), '+-');
      leg{l} = sprintf('Nb=%d',Nb(i)); l++;
    end
    axis([0 Fs/2 20 80]);
    hold off;
    h = legend(leg, "location", "northeast");
    legend("boxoff")

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit p-png", f);
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
      set (h, "fontsize", 16);
      xlabel('Freq (Hz)'); ylabel('Amplitude (dB)');
      print(sprintf("ratek2_nb_%s_%d",name,f),"-dpng","-S500,500");
    endif

  until (k == 'q')
  printf("\n");

endfunction
