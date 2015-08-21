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
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
%   $ cd ~/codec2-dev/octave
%   octave:14> newamp_fbf("../build_linux/src/hts1a",50)

function newamp_fbf(samname, f)
  
  newamp;

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);

  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
 
  plot_all_masks = 0;
  k = ' ';
  do 
    figure(1);
    clf;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    size(s);
    plot(s);
    axis([1 length(s) -20000 20000]);

    figure(2);
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);

    % plotting

    plot((1:L)*Wo*4000/pi, AmdB,";Am;r");
    axis([1 4000 0 80]);
    hold on;
    plot((0:255)*4000/256, Sw(f,:),";Sw;");

    [maskdB Am_freqs_kHz] = mask_model(AmdB, Wo, L);
    plot(Am_freqs_kHz*1000, maskdB, 'g');

    % Analysis by synthesis ---------------------------------------

    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    [newmaskdB local_maxima] = make_newmask(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz);

    plot(local_maxima(:,2)*Wo*4000/pi, 70*ones(1,length(local_maxima)), 'r+');
    plot(mask_sample_freqs_kHz*1000, newmaskdB, 'm');

    % optionally plot all masking curves

    if plot_all_masks
      mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
      for m=1:L
        maskdB = schroeder(m*Wo*4/pi, mask_sample_freqs_kHz) + AmdB(m);
        plot(mask_sample_freqs_kHz*1000, maskdB, "k--");
      end
    end

    hold off;

    % interactive menu

    printf("\rframe: %d  menu: n-next  b-back  p-png  q-quit m-all", f);
    fflush(stdout);
    k = kbhit();
    if (k == 'n')
      f = f + 1;
    endif
    if (k == 'b')
      f = f - 1;
    endif
    if k == 'm'
      if plot_all_masks == 0
         plot_all_masks = 1;
      else
         plot_all_masks = 0;
      end
    end
  until (k == 'q')
  printf("\n");

endfunction

