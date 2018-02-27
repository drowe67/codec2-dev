% c2wideband_fbf.m
%
% Copyright David Rowe 2017
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

% Interactive Octave script to explore frame by frame operation of
% wideband Codc 2.

#{
     Make sure codec2-dev is compiled with the -DDUMP option - see README for
     instructions.

     Usage:
      ~/codec2-dev/build_linux/src$ ./c2sim ~/Desktop/c2_hd/speech_orig_16k.wav --Fs 16000 --dump speech

      $ cd ~/codec2-dev/octave
      $ octave
      octave:1> c2wideband_fbf("../build_linux/src/speech")
#}

function c2wideband_fbf(samname, f=70, varargin)
  more off;
  newamp;

  Fs = 16000; Fs2 = Fs/2; K = 30;

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
    Am_freqs_kHz = (1:L)*Wo*Fs2/(1000*pi);

    [rate_K_vec rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model(f,:), K, Fs, 'para');

    % plots ----------------------------------
  
    figure(2); clf; 
    plot((1:L)*Wo*Fs2/pi, AmdB,";rate L;g+-");
    axis([1 Fs2 -20 80]);
    hold on;
    stem(rate_K_sample_freqs_kHz*1000, rate_K_vec, "b");

    [model_ AmdB_] = resample_rate_L(model(f,:), rate_K_vec, rate_K_sample_freqs_kHz, Fs, 'para');
    AmdB_ = AmdB_(1:L);

    plot((1:L)*Wo*Fs2/pi, AmdB_,";AmdB;r+-");
    hold off;

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit [%d]", f);
    fflush(stdout);
    k = kbhit();

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


