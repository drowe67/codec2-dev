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

function [decmaskdB masker_freqs_kHz min_error mse_log1 mse_log2] = make_decmask_abys(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz, freq_quant, amp_quant)

    Nsamples = 3;

    % search range

    m_st = 1;
    m_en = L;

    target = maskdB;
    dec_samples = [];
    mse_log1 = mse_log2 = zeros(Nsamples, L);

    % This step quite important.  Set noise floor to avoid washing machine/tinkling
    % sounds with bg noise as chunks of spectrum come and go.  The masking model
    % basis functions are not gd at representing continuous spectra, especially at 
    % LF.  Fortunately rather gd at voiced speech, which is a much tougher job.

    decmaskdB = 10*ones(1,L);

    for sample=1:Nsamples

      target_log = target;

      % find best position for sample to minimise distance to target

      min_mse = 1E32;

      for m=m_st:m_en
        single_mask_m = schroeder(m*Wo*4/pi, mask_sample_freqs_kHz,1) + AmdB(m);
        candidate = max(decmaskdB, single_mask_m);
        error = target - candidate;
        mse_log1(sample,m) = mse = sum(abs(error)); % MSE in log domain
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
          min_error = error;
        end
      end

      % Choose best candidate ------------------------------------------------

      decmaskdB = min_candidate;
      dec_samples = [dec_samples; AmdB(min_mse_m) min_mse_m];
      masker_freqs_kHz = dec_samples(:,2)*Wo*4/pi;
      %printf("sample: %d min_mse_m: %d\n", sample, min_mse_m);
    end

    bits = [];

    % sort into increasing freq order

    masker_amps_dB = dec_samples(:,1);
    masker_freqs_kHz = dec_samples(:,2)*Wo*4/pi;
    [fsrt fsrt_ind] = sort(masker_freqs_kHz);
    masker_freqs_kHz = fsrt;
    masker_amps_dB = masker_amps_dB(fsrt_ind);

    % Differential Freq Quantisers - sounds acceptable

    if freq_quant
           
      % first freq quant to harmonic number m=1:8

      f0_kHz = Wo*4/pi;
      [masker_freqs_kHz(1) abits] = quantise((1:8)*f0_kHz, masker_freqs_kHz(1));
      bits = [bits abits];
     
      % then quantise differences

      for i=2:Nsamples
        targ = masker_freqs_kHz(i) - masker_freqs_kHz(i-1);
        [q_freq abits] = quantise(0.2:0.2:1.6, targ);
        bits = [bits abits];
        masker_freqs_kHz(i) = masker_freqs_kHz(i-1) + q_freq;
      end

       decmaskdB = determine_mask(masker_amps_dB,  masker_freqs_kHz, mask_sample_freqs_kHz);
    end

    % Amplitude quantisation by fitting a straight line -------------------------
    % amp_quant == 1: high rate, quantise deltas
    % amp_quant == 2: low rate, don't quantise deltas

    if amp_quant > 0

      % Fit straight line

      f = masker_freqs_kHz*1000;
      [gradient intercept] = linreg(f, masker_amps_dB, Nsamples);
      % use quantised gradient to take into account quantisation
      % errors in rest of quantisation

      gradient_ = quantise(-0.04:0.005:0.04, gradient);
      %gradient_ = gradient;
      printf("gradient; %f gradient_: %f\n", gradient, gradient_);

      % determine deltas, or errors in straight line fit

      masker_amps_dB_lin = f*gradient_ + intercept;
      masker_amps_dB_lin_delta = masker_amps_dB - masker_amps_dB_lin;

      % optional plots

      if 0
        figure(10)
        clf;
        plot(f, masker_amps_dB, 'r+', 'markersize', 10, 'linewidth', 2)

        fplt = 0:100:3900
        hold on;
        plot(fplt, fplt*gradient + intercept, 'b')
        fplt*gradient + intercept

        % plot lines for deltas

        for i=1:length(f)
          y1 = f(i)*gradient + intercept;
          y2 = masker_amps_dB(i);
          plot([f(i) f(i)], [y1 y2], 'markersize', 10, 'linewidth', 2)
        end
        hold off;
      end

      % quantise the deltas

      masker_amps_dB_lin_delta_ = zeros(Nsamples,1);
      if amp_quant == 1
        for i=1:Nsamples
          masker_amps_dB_lin_delta_(i) = quantise(-21:3:21, masker_amps_dB_lin_delta(i));
          printf("dlin: %f dlin_: %f\n", masker_amps_dB_lin_delta(i), masker_amps_dB_lin_delta_(i));
          % masker_amps_dB_lin_delta_(i) = masker_amps_dB_lin_delta(i);
        end
      end

      masker_amps_dB = f*gradient_ + masker_amps_dB_lin_delta_ + intercept;
      %decmaskdB = determine_mask(masker_amps_dB,  masker_freqs_kHz, mask_sample_freqs_kHz);
      decmaskdB = determine_mask(masker_amps_dB,  masker_freqs_kHz, mask_sample_freqs_kHz);
    end

    if 0
    printf("\n");
    for i=1:Nsamples
      printf("freq: %f amp: %f\n", masker_freqs_kHz(i), masker_amps_dB(i));
    end
    end
endfunction


% quantise input sample to nearest value in table, optionally return bianry code

function [quant_out bits] = quantise(levels, quant_in)

  % find closest quantiser level

  best_se = 1E32;
  for i=1:length(levels)
    se = (levels(i) - quant_in)^2;
    if se < best_se
      quant_out = levels(i);
      best_se = se;
      best_i = i;
    end
  end

  % convert index to binary bits

  numbits = ceil(log2(length(levels)));
  bits = zeros(1, numbits);
  for b=1:numbits
    bits(b) = bitand(best_i-1,2^(numbits-b)) != 0;
  end

endfunction


function masker_freqs_kHz = unquantise_freqs(bits, Wo)
  for i=1:4
    st = (i-1)*3+1; en=i*3; 
    index(i) = bits(st:en) * [4 2 1]' + 1;
  end
  f0_kHz = Wo*4/pi;

  masker_freqs_kHz(1) = index(1)*f0_kHz;
     
  % then unquantise differences

  q_freqs = 0.2:0.2:1.6;

  for i=2:4
    masker_freqs_kHz(i) = masker_freqs_kHz(i-1) + q_freqs(index(i));
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

function maskdB = schroeder(freq_tone_kHz, mask_sample_freqs_kHz, bark_model=1)
  f_kHz = mask_sample_freqs_kHz;
  f_Hz = f_kHz*1000;

  % Schroeder Spreading Function

  if bark_model == 0
    dz = bark(freq_tone_kHz*1000)-bark(f_Hz);
  end

  if bark_model == 1

    % Modification by DR: Piecewise linear model that makes bands
    % beneath 1.5kHz wider to match the width of F1 and
    % "fill in" the spectra better for UV sounds.

    x1 = 0.5; x2 = 2;
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
  end

  if bark_model == 2

    % constant bandwidth model, useful for bg noise and UV

    %dz = bark(freq_tone_kHz*1000) - bark(f_Hz);
    dz = 0.2*bark(freq_tone_kHz*1000-f_Hz);
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

function maskdB_ = decimate_frame_rate(model, decimate, f, frames, mask_sample_freqs_kHz)
    max_amp = 80;

    Wo = model(f,1);
    L = min([model(f,2) max_amp-1]);

    % determine frames that bracket the one we need to interp

    left_f = decimate*floor((f-1)/decimate)+1; 
    right_f = left_f + decimate;
    if right_f > frames
      right_f = left_f;
    end

    % determine fraction of each frame to use

    left_fraction  = 1 - mod((f-1),decimate)/decimate;
    right_fraction = 1 - left_fraction;

    printf("f: %d left_f: %d right_f: %d left_fraction: %f right_fraction: %f \n", f, left_f, right_f, left_fraction, right_fraction)

    % fit splines to left and right masks

    left_Wo = model(left_f,1);
    left_L = min([model(left_f,2) max_amp-1]);
    left_AmdB = 20*log10(model(left_f,3:(left_L+2)));
    left_mask_sample_freqs_kHz = (1:left_L)*left_Wo*4/pi;

    right_Wo = model(right_f,1);
    right_L = min([model(right_f,2) max_amp-1]);
    right_AmdB = 20*log10(model(right_f,3:(right_L+2)));
    right_mask_sample_freqs_kHz = (1:right_L)*right_Wo*4/pi;

    maskdB_left_pp = splinefit(left_mask_sample_freqs_kHz, left_AmdB, left_L);
    maskdB_right_pp = splinefit(right_mask_sample_freqs_kHz, right_AmdB, right_L);

    % determine mask for left and right frames, sampling at Wo for this frame

    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    maskdB_left = ppval(maskdB_left_pp, mask_sample_freqs_kHz);
    maskdB_right = ppval(maskdB_right_pp, mask_sample_freqs_kHz);

    maskdB_ = left_fraction*maskdB_left + right_fraction*maskdB_right;
endfunction


% plot some masking curves, used for working on masking filter changes

function plot_masking
  Fs = 8000;

  figure(1)
  mask_sample_freqs_kHz = 0.1:0.1:(Fs/1000)/2;
  maskdB_cb = schroeder(0.5, mask_sample_freqs_kHz, 1);
  plot(mask_sample_freqs_kHz, maskdB_cb);
  hold on;
  maskdB_res = schroeder(0.5, mask_sample_freqs_kHz, 2);
  plot(mask_sample_freqs_kHz, maskdB_res,'g');

  for f=0.5:0.5:3
    maskdB_cb = schroeder(f, mask_sample_freqs_kHz, 1);
    plot(mask_sample_freqs_kHz, maskdB_cb);
    maskdB_res = schroeder(f, mask_sample_freqs_kHz, 2);
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

  [phase s] = mag_to_phase(Gdbfk);

endfunction


% AbyS returns f & a, this function plots values so we can consider quantisation

function plot_f_a_stats(f,a)

  % freq pdfs

  [fsrt fsrt_ind] = sort(f,2);
  fsrt /= 1000;
  figure(1)
  for i=1:4
    subplot(2,2,i)
    hist(fsrt(:,i),50)
    printf("%d min: %d max: %d\n", i, min(fsrt(:,i)), max(fsrt(:,i)))
    an_axis = axis;
    axis([0 4 an_axis(3) an_axis(4)])
  end

  % freq diff pdfs

  figure(2)
  for i=1:4
    subplot(2,2,i)
    if i == 1
      hist(fsrt(:,i),50)
    else
      hist(fsrt(:,i) - fsrt(:,i-1),50)
    end
    an_axis = axis;
    axis([0 4 an_axis(3) an_axis(4)])
  end

  % amplitude PDFs

  l = length(a);
  for i=1:l
    asrt(i,:) = a(i, fsrt_ind(i,:));
  end
  
  figure(3)
  for i=1:4
    subplot(2,2,i)
    hist(asrt(:,i) - mean(asrt(:,:),2))
    an_axis = axis;
    axis([-40 40 an_axis(3) an_axis(4)])
  end
  
  % find straight line fit

  for i=1:l
    [gradient intercept] = linreg(1000*fsrt(i,:), asrt(i,:), 4);
    alinreg(i,:) = gradient*1000*fsrt(i,:) + intercept;
    alinregres(i,:) = asrt(i,:) - alinreg(i,:);
    m(i) = gradient; c(i) = intercept;
  end

  figure(4)
  for i=1:4
    subplot(2,2,i)
    hist(alinregres(:,i))
    an_axis = axis;
    axis([-40 40 an_axis(3) an_axis(4)])
  end
  
  figure(5)
  subplot(211)
  m = m(find(m>-0.05));
  m = m(find(m<0.03));
  hist(m,50)
  title('gradient');
  subplot(212)
  c = c(find(c>0));
  hist(c,50)
  title('y-int');

endfunction
