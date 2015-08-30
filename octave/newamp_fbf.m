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
  
  more off;
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
    clf;
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);

    % plotting

    plot((1:L)*Wo*4000/pi, AmdB,";Am;r");
    axis([1 4000 0 80]);
    hold on;
    plot((1:L)*Wo*4000/pi, AmdB,";Am;r+");
    plot((0:255)*4000/256, Sw(f,:),";Sw;");

    [maskdB Am_freqs_kHz] = mask_model(AmdB, Wo, L);
    plot(Am_freqs_kHz*1000, maskdB, 'g');

    % optionally show harmonics that are not masked

    not_masked_m = find(maskdB < AmdB);
    if 0
      plot(not_masked_m*Wo*4000/pi, 70*ones(1,length(not_masked_m)), 'bk+');
    end

    % optionally plot synthesised spectrum (early simple model)

    if 0
      AmdB_ = maskdB;
      AmdB_(not_masked_m) += 6;
      plot(Am_freqs_kHz*1000, AmdB_, 'g');
      plot(Am_freqs_kHz*1000, AmdB_, 'g+');
    end

    % estimate low rate samples

    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    [decmaskdB local_maxima error_log candidate_log target_log] = make_decmask_abys(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz);
    
    [nlm tmp] = size(local_maxima(:,2));
    nlm = min(nlm,4);
    tonef_kHz = local_maxima(1:nlm,2)*Wo*4/pi;
    toneamp_dB = local_maxima(1:nlm,1);
    plot(tonef_kHz*1000, 70*ones(1,nlm), 'bk+');
    plot(mask_sample_freqs_kHz*1000, decmaskdB, 'm');

    % fit a line to amplitudes

    %[m b] = linreg(tonef_kHz, toneamp_dB, nlm);
    %plot(tonef_kHz*1000, tonef_kHz*m + b, "bk");
    %plot(tonef_kHz*1000, 60 + toneamp_dB - (tonef_kHz*m + b), "r+");

    % optionally plot all masking curves

    if plot_all_masks
      mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
      for m=1:L
        maskdB = schroeder(m*Wo*4/pi, mask_sample_freqs_kHz) + AmdB(m);
        plot(mask_sample_freqs_kHz*1000, maskdB, "k--");
      end
    end

    hold off;

    figure(3)
    plot(target_log,'g')
    hold on;
    plot(candidate_log(3,:),'b');
    plot(candidate_log(5,:),'b');
    plot(error_log(3,:),'r');
    plot(error_log(5,:),'r');
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

