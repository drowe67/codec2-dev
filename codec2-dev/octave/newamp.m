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

% dealing with UV, BG noise. Prob is flat spectra.  When we fit a low
% freq masking model it's formant shaped rather than flat.  So we get
% these gaps in spectra that come and go - waterfall noises.  In
% particular at low frequencies.  Good news is they don't need to be
% quantised too finely.  This model has the disadvantage of not having
% variable bandwidths.
% when do waterfall noises appear?
% idea: we could add the ability to have wider bands
%        Add some sort of slope or floor
%        increase spacing of samples?  Like min spacing in bark dimension

% can we fit a different shape?

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


% Alternative way to come up with a decimated mask model, using
% analysis by synthesis to determine the best place to put samples.
% Ahh, takes me back to when I was a slip of a speech coder, playing
% with my first CELP codec!

function [decmaskdB dec_samples error_log candidate_log target_log] = make_decmask_abys(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz)

    % band pass filter: limit search to 250 to 3800 Hz

    %m_st = max(1,floor((pi*250/4000)/Wo));
    %m_en = floor((pi*3800/4000)/Wo);
    m_st = 1;
    m_en = L;

    target = maskdB;
    decmaskdB = zeros(1,L);
    dec_samples = [];
    error_log = [];
    candidate_log = [];

    for sample=1:4

      target_log = target;

      % find best position for sample to minimise distance to target

      min_mse = 1E32;
      for m=m_st:m_en
        single_mask_m = schroeder(m*Wo*4/pi, mask_sample_freqs_kHz) + AmdB(m);
        candidate = max(decmaskdB, single_mask_m);
        error = target - candidate;
        mse = sum(abs(error)); % MSE in log domain
        error_log = [error_log; error];
        candidate_log = [candidate_log; candidate];
        %printf("m: %d f: %f error: %f\n", m, m*Wo*4/pi, mse);

        too_close = 0;
        for x=1:sample-1
          if abs(dec_samples(x,2) - m) < 2
            too_close = 1;
          end
        end

        if (mse < min_mse) && (too_close == 0)
          min_mse = mse;
          min_mse_m = m;
          min_candidate = candidate;
        end
      end

      decmaskdB = min_candidate;
      dec_samples = [dec_samples; AmdB(min_mse_m) min_mse_m];
      %printf("sample: %d min_mse_m: %d\n", sample, min_mse_m);
    end

    if 0
      % add floor to mask - unsuccessful attempt at fixing tinkle noises
      min_AmdB = max(dec_samples(:,1));
      decmaskdB = max(decmaskdB, min_AmdB - 20);
    end

    % simulate quantisation of amplitudes by adding some noise

    if 0
      masker_amps_dB = dec_samples(:,1);
      masker_amps_dB += 3*(1 - 2*rand(4,1));
      masker_freqs_kHz = dec_samples(:,2)*Wo*4/pi;
      decmaskdB = determine_mask(masker_amps_dB,  masker_freqs_kHz, mask_sample_freqs_kHz);
    end
    
    % quantisation of amplitudes.  Determine and subtract mean.  Quantise difference
    % from mean to 4 levels (2 bits), at 6dB/level:
    %
    %           Level
    %    0       -9            
    %    1       -3
    %    2       +3
    %    3       +9

    if 0
      masker_amps_dB = dec_samples(:,1);
      masker_freqs_kHz = dec_samples(:,2)*Wo*4/pi;

      energy_dB = mean(masker_amps_dB);
      masker_amps_dB -= energy_dB;
      for i=1:4
        masker_amps_dB(i) = quantise([-9 -3 3 9], masker_amps_dB(i));
      end
      masker_amps_dB += energy_dB;
      decmaskdB = determine_mask(masker_amps_dB,  masker_freqs_kHz, mask_sample_freqs_kHz);
    end

    % fit straight line to amplitudes (sounds bad, still a bug somewhere)

    if 0
      [gradient intercept] = linreg(dec_samples(:,2), dec_samples(:,1), 4);
      masker_amps_dB = dec_samples(:,2)*gradient + intercept;
      masker_freqs_kHz = dec_samples(:,2)*Wo*4/pi;
      decmaskdB = determine_mask(masker_amps_dB,  masker_freqs_kHz, mask_sample_freqs_kHz);
    end
endfunction


% quantise input sample to nearest value in table 

function quant_out = quantise(levels, quant_in)
  best_se = 1E32;
  for i=1:length(levels)
    se = (levels(i) - quant_in)^2;
    if se < best_se
      quant_out = levels(i);
      best_se = se;
    end
  end
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


% Sample mask as model for Am

function [maskdB Am_freqs_kHz] = mask_model(AmdB, Wo, L)

    Am_freqs_kHz = (1:L)*Wo*4/pi;
    maskdB = determine_mask(AmdB, Am_freqs_kHz, Am_freqs_kHz);
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

    x1 = 0.5; x2 = 1.5;
    y1 = 0.5; y2 = 1;
    grad  = (y2 - y1)/(x2 - x1);
    y_int = y1 - grad*x1;

    if freq_tone_kHz <= x1
      y = y1;
    end
    if (freq_tone_kHz > x1) && (freq_tone_kHz < x2)
      y = grad*freq_tone_kHz + y_int;
    end
    if freq_tone_kHz >= x2
      y = y2;
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


% -12dB/octave mask to model speech articulation

function maskdB = resonator(freq_tone_kHz, mask_sample_freqs_kHz)
  maskdB = zeros(1, length(mask_sample_freqs_kHz));
  for m=1:length(mask_sample_freqs_kHz)
    maskdB(m) = -24*abs(log2(freq_tone_kHz/mask_sample_freqs_kHz(m)));
    %printf("m: %d ft: %f fm: %f ft/fm: %f maskdB: %f\n", m, freq_tone_kHz, mask_sample_freqs_kHz(m), freq_tone_kHz/mask_sample_freqs_kHz(m), maskdB(m));
  end
endfunction


% sampling the mask in one frame using an arbitrary set of samplng frequencies

function maskdB = resample_mask(model, f, mask_sample_freqs_kHz)
  max_amp = 80;

  Wo = model(f,1);
  L = min([model(f,2) max_amp-1]);
  Am = model(f,3:(L+2));
  AmdB = 20*log10(Am);
  masker_freqs_kHz = (1:L)*Wo*4/pi;
  maskdB = determine_mask(AmdB, masker_freqs_kHz, mask_sample_freqs_kHz);
endfunction


% decimate frame rate of mask, use linear interpolation in the log domain 

function maskdB_ = decimate_frame_rate(maskdB, model, decimate, f, frames, mask_sample_freqs_kHz)
    max_amp = 80;

    Wo = model(f,1);
    L = min([model(f,2) max_amp-1]);

    if (mod(f-1,decimate) && (f < (frames-decimate)))

      % determine frames that bracket the one we need to interp

      left_f = decimate*floor(f/decimate)+1; 
      right_f = decimate*ceil(f/decimate)+1;

      % determine fraction of each frame to use

      right_fraction  = mod(f-1,decimate)/decimate;
      left_fraction = 1 - right_fraction;

      printf("\nf: %d left_f: %d right_f: %d left_fraction: %f right_fraction: %f \n",f,left_f,right_f,left_fraction,right_fraction)

      % determine mask for left and right frames, sampling at Wo for this frame

      mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
      maskdB_left = resample_mask(model, left_f, mask_sample_freqs_kHz);
      maskdB_right = resample_mask(model, right_f, mask_sample_freqs_kHz);

      maskdB_ = left_fraction*maskdB_left + right_fraction*maskdB_right;
    else
      maskdB_ = maskdB;
      printf("\n");
    end
endfunction


% plot some masking curves, used for working on masking filter changes

function plot_masking
  Fs = 8000;

  figure(1)
  mask_sample_freqs_kHz = 0.1:0.1:(Fs/1000)/2;
  maskdB_cb = schroeder(0.5, mask_sample_freqs_kHz, 1);
  plot(mask_sample_freqs_kHz, maskdB_cb);
  hold on;
  maskdB_res = resonator(0.5, mask_sample_freqs_kHz);
  plot(mask_sample_freqs_kHz, maskdB_res,'g');

  for f=0.5:0.5:3
    maskdB_cb = schroeder(f, mask_sample_freqs_kHz, 1);
    plot(mask_sample_freqs_kHz, maskdB_cb);
    maskdB_res = resonator(f, mask_sample_freqs_kHz);
    plot(mask_sample_freqs_kHz, maskdB_res,'g');
  end
  hold off;
  axis([0.1 4 -30 0])
  grid

  figure(2)
  clf;
  w = pi/4; beta = 0.9;
  X = freqz(1,[1 -2*beta*cos(w) beta*beta],4000);
  plot(10*log10(abs(X)))
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


% Determine a phase spectra from a magnitude spectra
% from http://www.dsprelated.com/showcode/20.php
% Haven't _quite_ figured out how this works but have to start somewhere ....
%
% TODO: we may be able to sample at a lower rate, like mWo
%       but start with something that works

function [phase Gdbfk s Aw] = determine_phase(model, f, ak)
  Nfft    = 512;  % FFT size to use 
  Fs      = 8000;
  max_amp = 80;
  L       = min([model(f,2) max_amp-1]);
  Wo      = model(f,1);

  mask_sample_freqs_kHz = (Fs/1000)*[0:Nfft/2]/Nfft;           % fft frequency grid (nonneg freqs)
  Gdbfk = resample_mask(model, f, mask_sample_freqs_kHz);

  % optional input of aks for testing

  if nargin == 3
    Aw = 1 ./ fft(ak,Nfft);
    Gdbfk = 20*log10(abs(Aw(1:Nfft/2+1)));
  end

  Ns = length(Gdbfk); if Ns~=Nfft/2+1, error("confusion"); end
  Sdb = [Gdbfk,Gdbfk(Ns-1:-1:2)]; % install negative-frequencies

  S = 10 .^ (Sdb/20); % convert to linear magnitude
  s = ifft(S); % desired impulse response
  s = real(s); % any imaginary part is quantization noise
  tlerr = 100*norm(s(round(0.9*Ns:1.1*Ns)))/norm(s);
  disp(sprintf(['  Time-limitedness check: Outer 20%% of impulse ' ...
               'response is %0.2f %% of total rms'],tlerr));
  % = 0.02 percent

  if tlerr>1.0 % arbitrarily set 1% as the upper limit allowed
    disp('  Increase Nfft and/or smooth Sdb\n');
  end

  c = ifft(Sdb); % compute real cepstrum from log magnitude spectrum

  % Check aliasing of cepstrum (in theory there is always some):

  caliaserr = 100*norm(c(round(Ns*0.9:Ns*1.1)))/norm(c);
  disp(sprintf(['  Cepstral time-aliasing check: Outer 20%% of ' ...
               'cepstrum holds %0.2f %% of total rms\n'],caliaserr));

  if caliaserr>1.0 % arbitrary limit
    disp('  Increase Nfft and/or smooth Sdb to shorten cepstrum\n');
  end

  % Fold cepstrum to reflect non-min-phase zeros inside unit circle:

  cf = [c(1), c(2:Ns-1)+c(Nfft:-1:Ns+1), c(Ns), zeros(1,Nfft-Ns)];
  Cf = fft(cf); % = dB_magnitude + j * minimum_phase

  phase = imag(Cf)/(20/log(10));
endfunction

