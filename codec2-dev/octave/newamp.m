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
%       [ ] consider a floor in mask to interpolate missing bits
%       [ ] some sort of filter on min amp/diff from previous peak, e.g. vk5qi:161, 172
%           + min spacing in bark?  Make larger at higher freq?
%           + or some other measure, like making sure choice minimises MSE
%   [ ] way to output at various processing steps, like PF initial mask, pre PF
%   [ ] BPF at all?
%   [ ] phase model?  Fit LPC, just swing phase at peaks? Try no phase tweaks

1;

% Create a "decimated mask" model using just a few samples of
% critical band filter centre frequencies.  For voiced speech,
% we fit the amplitude of these samples to a straight line.
% TODO: 
%  [ ] track down bg noises on vk5qi and kristoff
%  [ ] return data for plotting, like slope m
%  [ ] quantise m

function [decmaskdB local_maxima_sort] = make_decmask(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz)

    % band pass filter: limit search to 250 to 3800 Hz

    m_st = max(1,floor((pi*250/4000)/Wo));
    m_en = floor((pi*3800/4000)/Wo);

    % We start off by assuming that the local maxima in the masking
    % curve are the centres of the samples we want to keep.

    local_maxima = [];
    if maskdB(m_st) > maskdB(m_st+1)
      local_maxima = [local_maxima; AmdB(m_st) m_st];
    end
    for m=m_st+1:m_en-1
      if (maskdB(m-1) < maskdB(m)) && (maskdB(m) > maskdB(m+1))
        local_maxima = [local_maxima; AmdB(m) m];
      end
    end
    [nlm tmp] = size(local_maxima);

    % Occasionally there are no local maxima so pop one in

    if nlm == 0
      local_maxima = [AmdB(m_st) m_st];
      nlm = 1;
    end
    
    % fit a straight line to the amplitudes of our candidate samples,
    % this will help us later when we code and transmit the amplitude 
    % of each sample

    if nlm > 1
      [m b] = linreg(local_maxima(:,2), local_maxima(:,1), nlm);
      local_maxima_fit = local_maxima(:,2)*m + b;
    else
      local_maxima_fit = local_maxima(1,1);
    end

    % Remove any outliers to the straight line fit: Sometimes local
    % maxima appear in an antiformant regions, say if F1 and F2 are a
    % long way apart.  After a straight line fit the anti-format
    % amplitide sample will be way off the straight line, which will
    % cause a spike of spectral energy right where we don't want it -
    % in the middle of an antiformat.  So lets test the fit of each
    % sample, and only include those that work well with the straight
    % line fit.  For voiced frames, m < 0.  For UV frames, we don't
    % care about the straight line fit as unvoiced speech is all over
    % the place in amplitude anyway.
    
    local_maxima2 = [];
    for i=1:nlm
      if (local_maxima_fit(i) - local_maxima(i,1) < 6) || (m > 0)
        local_maxima2 = [local_maxima2; local_maxima(i,1) local_maxima(i,2)];
      end
    end

    % now sort and keep the top 4 samples
    
    local_maxima_sort = flipud(sortrows(local_maxima2,1));
    [nlm tmp] = size(local_maxima_sort);
    nlm = min(nlm,4);
    local_maxima_sort = local_maxima_sort(1:nlm,:);

    % fit straight line again, this time with outliers removed

    [m b] = linreg(local_maxima_sort(:,2), local_maxima_sort(:,1), nlm);
    masker_amps_dB = local_maxima_sort(:,2)*m + b;
    masker_freqs_kHz = local_maxima_sort(:,2)*Wo*4/pi;
    %masker_amps_dB = local_maxima_sort(:,1);
    %masker_freqs_kHz = local_maxima_sort(:,2)*Wo*4/pi;

    % and construct new, decimated mask using our small set of
    % samples, with amplitudes fitted to a linear line

    decmaskdB = determine_mask(masker_amps_dB,  masker_freqs_kHz, mask_sample_freqs_kHz);
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


% produce a scatter diagram of amplitudes

function amp_scatter(samname)
  
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);
  max_amp = 80;

  Am_out_name = sprintf("%s_am.out", samname);
  freqs = [];
  ampsdB = [];

  for f=1:frames
    
    L = min([model(f,2) max_amp-1]);
    Wo = model(f,1);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);

    maskdB = mask_model(AmdB, Wo, L);
    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    [newmaskdB local_maxima] = make_newmask(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz);

    [nlm tmp] = size(local_maxima);
    freqs = [freqs (local_maxima(1:min(4,nlm),2)*Wo*4000/pi)'];
    an_ampsdB = local_maxima(1:min(4,nlm),1)';
    ampsdB = [ampsdB an_ampsdB-mean(an_ampsdB)];
  end

  figure(1)
  plot(freqs, ampsdB,'+');
  figure(2)
  subplot(211)
  hist(freqs,20)
  subplot(212)
  hist(ampsdB,20)
endfunction

%amp_scatter("../build_linux/src/all")
