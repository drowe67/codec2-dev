% newamp.m
%
% Copyright David Rowe 2015
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Library of Octave functions to explore new ideas in amplitude
% (spectral envelope) modelling.  See newamp_fby (frame by frame
% analysis) and newamp_batch (batch processing for listening tests)
%
% we don't care about
%  + spectral tilt, in can vary on input and our model shouldnt care.
%    We can vary it on output and the listener won't care 
%  + absolute energy of entire signal
%  + harmonics beneath the masking curve
% we do care about:
%  + clearly defined formant formation
%  + some relative amplitude info, like dominance of HF for UV sounds
%
% TODO:
%   [ ] waterfall sounds
%       [X] tweak CB bandwidths at LF to be wider
%   [ ] way to output at various processing steps, like PF initial mask, pre PF
%   [ ] BPF at all?
%   [ ] phase model?  Fit LPC, just swing phase at peaks? Try no phase tweaks

1;

% model mask using just a few samples of critical band filter centre frequencies

function [newmaskdB local_maxima_sort] = make_newmask(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz)

    % find local maxima and sort in descending order of magnitude

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

    % construct new mask from a small subset of sorted local maxima,
    % this is the highly compressed model of the amplitude envelope
    % that we send to the decoder

    masker_amps_dB = local_maxima_sort(1:min(4,nlm),1);
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


%
% Masking functions from http://www.perceptualentropy.com/coder.html#C
% Thanks Jon Boley!
%

% Calculate the Schroeder masking spectrum for a given frequency and SPL

function maskdB = schroeder(freq_tone_kHz, mask_sample_freqs_kHz, modified_bark_en=1)
  f_kHz = mask_sample_freqs_kHz;
  f_Hz = f_kHz*1000;

  % Schroeder Spreading Function

  if modified_bark_en == 1

    % Modification by DR: Piecewise linear model that makes bands
    % beneath 1.5kHz wider to match the width of F1 and
    % "fill in" the spectra better for UV sounds.

    if freq_tone_kHz <= 0.5
      y = 0.5;
    end
    if (freq_tone_kHz > 0.5) && (freq_tone_kHz < 1.5)
      y = 0.5*freq_tone_kHz + 0.25;
    end
    if freq_tone_kHz >= 1.5
      y = 1;
    end
    dz = y*(bark(freq_tone_kHz*1000) - bark(f_Hz));
  else
    dz = bark(freq_tone_kHz*1000)-bark(f_Hz);
  end

  maskdB = 15.81 + 7.5*(dz+0.474) - 17.5*sqrt(1 + (dz+0.474).^2);
endfunction


% Converts frequency to bark scale
% Frequency should be specified in Hertz

function b=bark(f)
  b = 13*atan(0.76*f/1000) + 3.5*atan((f/7500).^2); 
endfunction


% plot some masking curves, used for working on masking filter changes

function plot_masking
  Fs = 8000;

  figure(1)
  mask_sample_freqs_kHz = 0.1:0.1:(Fs/1000)/2;
  maskdB = schroeder(0.5, mask_sample_freqs_kHz, 0);
  plot(mask_sample_freqs_kHz, maskdB);
  hold on;
  maskdB = schroeder(0.5, mask_sample_freqs_kHz, 1);
  plot(mask_sample_freqs_kHz, maskdB,'g');

  for f=1:0.5:3
    maskdB = schroeder(f, mask_sample_freqs_kHz, 0);
    plot(mask_sample_freqs_kHz, maskdB);
    maskdB = schroeder(f, mask_sample_freqs_kHz, 1);
    plot(mask_sample_freqs_kHz, maskdB,'g');
  end
  hold off;
  axis([0.1 4 -30 0])

  figure(2)
  plot(mask_sample_freqs_kHz, bark(mask_sample_freqs_kHz*1000))
  hold on;
  plot(mask_sample_freqs_kHz, modified_bark(mask_sample_freqs_kHz*1000),'g')
  hold off;
  grid
endfunction

