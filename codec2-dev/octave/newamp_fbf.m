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
  
  more off;
  newamp;
  phase_stuff = 0;
  plot_not_masked = 0;
  plot_spectrum = 1;
  freq_quant = 0;
  amp_quant = 0;

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);

  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);

  if phase_stuff
    ak_name = strcat(samname,"_ak.txt");
    ak = load(ak_name);
  end

  pp_bw = gen_pp_bw;

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
    %[h w] = freqz([1 -1],1,(1:L)*Wo); % pre-emphasise to reduce dynamic range
    %Am = Am .* abs(h);
    AmdB = 20*log10(Am);

    % plotting

    axis([1 4000 -20 80]);
    hold on;
    if plot_spectrum
      plot((1:L)*Wo*4000/pi, AmdB,";Am;r");
      plot((1:L)*Wo*4000/pi, AmdB,"r+");
      %plot((0:255)*4000/256, Sw(f,:),";Sw;");
    end

    [maskdB Am_freqs_kHz] = mask_model(AmdB, Wo, L);

    %plot(Am_freqs_kHz*1000, maskdB, 'g');
    %plot(Am_freqs_kHz2*1000, maskdB2, 'r');

    % optionally plot synthesised spectrum (early simple model)

    if 0
      AmdB_ = maskdB;
      plot(Am_freqs_kHz*1000, AmdB_, 'g');
      plot(Am_freqs_kHz*1000, AmdB_, 'g+');
    end

    % decimate in frequency

    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    [decmaskdB masker_freqs_kHz min_error mse_log1 mse_log2] = make_decmask_abys(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz, freq_quant, amp_quant, pp_bw);

    % find turning points - prototype for finding PF freqs when we decimate in time
    
    d = decmaskdB(2:L) - decmaskdB(1:L-1);
    tp = [];
    for m=1:L-1
      if (d(m) > 0) && (d(m+1) < 0)
        tp = [tp m+1];
      end
    end

    figure(2)
    
    tonef_kHz = masker_freqs_kHz;
    nlm = length(tonef_kHz);

    plot(tonef_kHz*1000, zeros(1,nlm), ';AbyS Mask Freqs;bk+');
    plot(mask_sample_freqs_kHz*1000, decmaskdB, ';AbyS Mask;m');
    plot(tp*Wo*4000/pi, 5*ones(1,length(tp)), ';PF Freqs;r+');
    legend('boxoff');
    %plot(mask_sample_freqs_kHz*1000, min_error);
    hold off;

    figure(3)
    clf
    plot((1:L)*Wo*4000/pi, mse_log1');
    axis([0 4000 0 max(mse_log1(1,:))])
    title('Error as a function of position for each AbyS stage');

    % decimated in time

    %maskdB = decimate_frame_rate(maskdB, model, 4, f, frames, mask_sample_freqs_kHz);
    %plot(mask_sample_freqs_kHz*1000, maskdB, 'k');

    % optionally plot all masking curves

    if plot_all_masks
      mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
      for m=1:L
        maskdB = schroeder(m*Wo*4/pi, mask_sample_freqs_kHz) + AmdB(m);
        plot(mask_sample_freqs_kHz*1000, maskdB, "k--");
      end
    end

    hold off;

    if phase_stuff

      [phase Sdb s Aw] = determine_phase(model, f, ak(f,:));
      figure(3)
      subplot(211)
      plot(Sdb)
      title('Mag (dB)');
      subplot(212)
      plot(phase(1:256))
      hold on;
      plot(angle(Aw(1:256))+0.5,'g')
      hold off;
      title('Phase (rads)');
      figure(4)
      plot(s)
    end

    % interactive menu

    printf("\rframe: %d  menu: n-next  b-back a-Am ", f);
    if freq_quant
      printf("F");
    else
      printf("f");
    end
    printf("-freq_quant ");
    if amp_quant
      printf("M");
    else
      printf("m");
    end
    printf("-amp_quant q-quit");

    fflush(stdout);
    k = kbhit();
    if (k == 'n')
      f = f + 1;
    endif
    if (k == 'b')
      f = f - 1;
    endif
    if k == 'm'
      if amp_quant == 0
         amp_quant = 1;
      else
         amp_quant = 0;
      end
    end
    if k == 'f'
      if freq_quant == 0
        freq_quant = 1;
      else
        freq_quant = 0;
      end
    end
    if k == 'a'
      if plot_spectrum == 0
         plot_spectrum = 1;
      else
         plot_spectrum = 0;
      end
    end
  until (k == 'q')
  printf("\n");

endfunction

