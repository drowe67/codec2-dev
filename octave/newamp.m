% Copyright David Rowe 2015
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Octave script to explore new ideas in amplitude modelling.
%
% we don't care about
%  + spectral tilt, in can vary on input and our quantiser shouldnt care.
%    We can vary it on output and the listener won't care 
%  + absolute energy of entire signal
%  + harmonics beneath the masking curve
% we do care about:
%  + clearly defined formant formation
%  + some relative amplitude info, like dominance of HF for UV sounds
%
% TODO:
%   [ ] waterfall sounds
%       [ ] tweak CB bandwidths at LF to be wider
%   [ ] way to output at various processing steps, like PF initial mask, pre PF
%   [ ] BPF at all?
%   [ ] phase model?  Fit LPC, just swing phase at peaks? Try no phase tweaks

1;

% 
% frame by frame analysis

function newamp_frame(samname, f)
  
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


% model mask using just a few samples of critical band filters

function [newmaskdB local_maxima_sort] = make_newmask(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz)
    local_maxima = [];
    if maskdB(1) > maskdB(2)
      local_maxima = [local_maxima; AmdB(1) 1];
    end
    for m=2:L-1
      if (maskdB(m-1) < maskdB(m)) && (maskdB(m) > maskdB(m+1))
        local_maxima = [local_maxima; AmdB(m) m];
      end
    end
    [nlm tmp] = size(local_maxima);
    if nlm == 0
      local_maxima = [AmdB(1) 1];
      nlm = 1;
    end
    local_maxima_sort = flipud(sortrows(local_maxima,1));

    % construct new mask from sorted local maxima

    masker_amps_dB = local_maxima_sort(1:min(4,nlm),1);
    % masker_amps_dB = mean(masker_amps_dB)*ones(1,min(4,nlm));
    masker_freqs_kHz = local_maxima_sort(1:min(4,nlm),2)*Wo*4/pi;
    newmaskdB = determine_mask(masker_amps_dB, masker_freqs_kHz, mask_sample_freqs_kHz);
endfunction


% determine cumulative mask, using amplitude of each harmonic.  Mask is
% sampled across L points in the linear domain

function maskdB = determine_mask(masker_amps_dB, masker_freqs_kHz, mask_sample_freqs_kHz)

    % calculate and plot masking curve

    maskdB = zeros(1,length(mask_sample_freqs_kHz));
    for m=1:length(masker_freqs_kHz)
      maskdB = max(maskdB, schroeder(masker_freqs_kHz(m), mask_sample_freqs_kHz) + masker_amps_dB(m)); 
    end
end


% Sample mask as model for Am, tweaked a bit to enhance
% the formants:antiformant radtio, like the LPC postfilter

function [maskdB_pf Am_freqs_kHz peaks_kHz] = mask_model(AmdB, Wo, L)

    Am_freqs_kHz = (1:L)*Wo*4/pi;
    maskdB = determine_mask(AmdB, Am_freqs_kHz, Am_freqs_kHz);
    non_masked_m = find(maskdB < AmdB);
    peaks_kHz = non_masked_m*Wo*4/pi;
    maskdB_pf = maskdB;
    %maskdB_pf = zeros(1,L);
    %maskdB_pf(non_masked_m) = maskdB(non_masked_m) + 6;

endfunction


% process a whole file and write results

function newamp_batch(samname)
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);
  max_amp = 80;

  Am_out_name = sprintf("%s_am.out", samname);
  fam = fopen(Am_out_name,"wb"); 

  for f=1:frames
    
    L = min([model(f,2) max_amp-1]);
    Wo = model(f,1);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);

    maskdB = mask_model(AmdB, Wo, L);
    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    [newmaskdB local_maxima] = make_newmask(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz);

    [nlm tmp] = size(local_maxima);
    non_masked_m = local_maxima(1:min(4,nlm),2);
    maskdB_pf = newmaskdB;
    maskdB_pf(non_masked_m) = maskdB_pf(non_masked_m) + 6;
    
    Am_ = zeros(1,max_amp);
    Am_(2:L) = 10 .^ (maskdB_pf(1:L-1)/20);
    
    % save to file
    fwrite(fam, Am_, "float32");
  end

  fclose(fam);

endfunction

%
% Masking functions from http://www.perceptualentropy.com/coder.html#C
% Thanks Jon Boley!
%

% Calculate the Schroeder masking spectrum for a given frequency and SPL

function maskdB = schroeder(freq_tone_kHz, mask_sample_freqs_kHz)
  f_kHz = mask_sample_freqs_kHz;
  A = 3.64*(f_kHz).^(-0.8) - 6.5*exp(-0.6*(f_kHz - 3.3).^2) + (10^(-3))*(f_kHz).^4;
  f_Hz = f_kHz*1000;

  % Schroeder Spreading Function
  dz = bark(freq_tone_kHz*1000)-bark(f_Hz);
  maskdB = 15.81 + 7.5*(dz+0.474) - 17.5*sqrt(1 + (dz+0.474).^2);
endfunction


% Converts frequency to bark scale
% Frequency should be specified in Hertz

function b=bark(f)
  b = 13*atan(0.76*f/1000) + 3.5*atan((f/7500).^2); 
endfunction


% return a sampling grid of frequencies in Hz given B equally spaced
% points on the bark scale

function f_Hz = bark_freq_samples(B)
   Fs2 = 4000;
   bark_lut = bark(1:Fs2);
   bark_grid_step = 1/B;
   bark_grid = (bark_grid_step:bark_grid_step:1)*bark(Fs2);
   f_Hz = interp1(bark_lut, 1:Fs2, bark_grid);
endfunction


% plot some masking curves

function plot_masking
  Fs = 8000;
  mask_sample_freqs_kHz = 0.1:0.1:(Fs/1000)/2;
  maskdB = schroeder(1, mask_sample_freqs_kHz)
  figure(1)
  plot(mask_sample_freqs_kHz, maskdB);
endfunction

% Choose one of these to run

more off;
%newamp_frame("../build_linux/src/hts1a", 121);
newamp_batch("../build_linux/src/hts1a");
%plot_masking
