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
  decimate_in_freq = 0;

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

  % pp_bw = gen_pp_bw;

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

    axis([1 4000 -20 80]);
    hold on;
    if plot_spectrum
      plot((1:L)*Wo*4000/pi, AmdB,";Am;r");
      plot((1:L)*Wo*4000/pi, AmdB,"r+");
      %plot((0:255)*4000/256, Sw(f,:),";Sw;");
    end

    [maskdB Am_freqs_kHz] = mask_model(AmdB, Wo, L);

    plot(Am_freqs_kHz*1000, maskdB, 'g');

    % Lets try to come up with a smoothed model.  Replace points from 3500 Hz on with
    % those than match low end of spectrum.  This will make it more cyclical and make
    % DFT happier, less high freq energy.  Need a better way to describe that.

    anchor = floor(7*L/8);
    xpts = [ anchor-1 anchor L+1 L+2];
    ypts = [ maskdB(anchor-1) maskdB(anchor) maskdB(1) maskdB(2)];
    mask_pp = splinefit(xpts, ypts, 1);
    maskdB_smooth = [maskdB(1:anchor) ppval(mask_pp, anchor+1:L)];

    plot(Am_freqs_kHz*1000, maskdB_smooth, 'b');

    if decimate_in_freq
      % decimate in frequency

      mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
      [decmaskdB masker_freqs_kHz masker_amps_dB min_error mse_log1] = make_decmask_mel(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz, freq_quant, amp_quant);

    % find turning points - prototype for finding PF freqs when we decimate in time
    
    if 0
      d = decmaskdB(2:L) - decmaskdB(1:L-1);
      size(d)
      tp = [];
      for m=1:L-1
        if (d(m) > 0) && (d(m+1) < 0)
          tp = [tp m+1];
        end
      end
    end

      figure(2)
    
      tonef_kHz = masker_freqs_kHz;
      nlm = length(tonef_kHz);

      plot(tonef_kHz*1000, zeros(1,nlm), ';AbyS Mask Freqs;bk+');
      if decimate_in_freq
        plot(mask_sample_freqs_kHz*1000, decmaskdB, ';AbyS Mask;m');
      end
      %plot(tp*Wo*4000/pi, 5*ones(1,length(tp)), ';PF Freqs;r+');
      legend('boxoff');
      %plot(mask_sample_freqs_kHz*1000, min_error);

    end

    hold off;

    % lets get a feel for the "spectrum" of the smoothed spectral envelope
    % this will give us a feel for how hard it is to code, ideally we would like
    % just a few coefficents to be non-zero

    figure(3)
    clf

    D = abs(fft(maskdB));
    en = floor(length(D)/2+1);
    stem(D(2:en),'g')
    hold on;
    D_smooth = abs(fft(maskdB_smooth));
    stem(D_smooth(2:en),'b')
    hold off;

    % let plot the cumulative amount of energy in each DFT

    figure(4)
    clf
    plot(cumsum(D(2:en)/sum(D(2:en))),'g');
    hold on;
    plot(cumsum(D_smooth(2:en)/sum(D_smooth(2:en))),'b');
    hold off;
    axis([1 L 0 1])

    % try truncating DFT coeffs and see what that looks like

    D = fft(maskdB_smooth);
    keep = 7;
    D(keep+1:L-keep) = 0;
    d = ifft(D);

    figure(2)
    hold on;
    plot(Am_freqs_kHz*1000, real(d), 'c')
    hold off;

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
         amp_quant = 3;
      else
         amp_quant = 0;
      end
    end
    if k == 'f'
      if freq_quant == 0
        freq_quant = 3;
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

