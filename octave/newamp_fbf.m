% newamp_fbf.m
%
% Copyright David Rowe 2015
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
%   octave:14> newamp_fbf("../build_linux/src/hts1a",50)



function newamp_fbf(samname, f=10)
  newamp;
  more off;
  plot_spectrum = 1;
  vq_en = 1;

  % load up text files dumped from c2sim ---------------------------------------

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);
  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);

  if vq_en
    load vq;
  end

  % Keyboard loop --------------------------------------------------------------

  k = ' ';
  do 
    figure(1);
    clf;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    size(s);
    plot(s);
    axis([1 length(s) -20000 20000]);

    figure(2);
    clf;
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);

    axis([1 4000 -20 80]);
    hold on;
    if plot_spectrum
      plot((1:L)*Wo*4000/pi, AmdB,";Am;r");
      plot((1:L)*Wo*4000/pi, AmdB,"r+");
    end

    [maskdB Am_freqs_kHz] = mask_model(AmdB, Wo, L);
    [tmp1 tmp2 D] = decimate_in_freq(maskdB, 0);
    if vq_en
      [maskdB_ maskdB_cyclic D_cyclic dk_] = decimate_in_freq(maskdB, 1, 7, vq);
    else
      [maskdB_ maskdB_cyclic D_cyclic dk_] = decimate_in_freq(maskdB, 1);
    end

    plot(Am_freqs_kHz*1000, maskdB, ';mask;g');
    plot(Am_freqs_kHz*1000, maskdB_cyclic, ';mask cyclic;b');
    plot(Am_freqs_kHz*1000, maskdB_, ';mask trunated;c');

    % lets get a feel for the "spectrum" of the smoothed spectral envelope
    % this will give us a feel for how hard it is to code, ideally we would like
    % just a few coefficents to be non-zero

    figure(3)
    clf

    en = L/2+1;
    stem(D(2:en),'g')
    hold on;
    stem(D_cyclic(2:en),'b')
    hold off;

    % let plot the cumulative amount of energy in each DFT

    figure(4)
    clf
    plot(cumsum(D(2:en)/sum(D(2:en))),';cumsum;g');
    hold on;
    plot(cumsum(D_cyclic(2:en)/sum(D_cyclic(2:en))),';cumsum cyclic;b');
    hold off;
    axis([1 L 0 1])

    figure(5)
    clf
    stem(dk_)

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back q-quit", f);
    fflush(stdout);
    k = kbhit();

    if (k == 'n')
      f = f + 1;
    endif
    if (k == 'b')
      f = f - 1;
    endif
  until (k == 'q')
  printf("\n");

endfunction


