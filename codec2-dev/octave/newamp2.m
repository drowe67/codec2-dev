% newamp2.m
%
% Copyright David Rowe 2018
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Library of Octave functions to explore new ideas in amplitude
% (spectral envelope) modelling.  See newamp2_fby (frame by frame
% analysis) and newamp2_batch (batch processing for listening tests)
%

1;

% --------------------------------------------------------------------------------
% Functions used by rate K mel work
% --------------------------------------------------------------------------------

function y = lanczos2(x)
  y = sinc(x).*sinc(x/2);
endfunction

function y = interp_lanczos(xp, yp, xp_max, x)
  
  y = zeros(1,length(x));
  k = 1;
  for i=1:length(x)
    % find closest sample in xp just greater than xi
    xi = x(i);
    while (xp(k) <= xi) && (k < length(xp)-2)
      k++;
    end
    
    % we'd like to use k-2 .. k+2, but we need to stay inside limits of xp

    k_st = k - 2; k_st = max(1,k_st);
    k_en = k + 2; k_en = min(length(xp),k_en);
    % printf("i: %d xi: %f k: %d k_st: %d k_en: %d\n", i, xi, k, k_st, k_en);
    
    % map frequencies to x in -2 ... + 2

    delta = xp(2) - xp(1);
    xl = (xp(k_st:k_en) - xi)/delta;

    y(i) = lanczos2(xl) * yp(k_st:k_en)';
    
  end
  
endfunction


% General 2nd order parabolic interpolator.  Used splines orginally,
% but this is much simpler and we don't need much accuracy.  Given two
% vectors of points xp and yp, find interpolated values y at points x
%

% If a point in x is less than the smallest point in xp, we linearly
% interpolate down to (0,0).  If a point in x is greater than the
% greatest value in xp, we linearly interpolate down to (xp_max, 0)

function y = interp_para(xp, yp, xp_max, x)
  assert( (length(xp) >=3) && (length(yp) >= 3) );

  y = zeros(1,length(x));
  k = 1;
  for i=1:length(x)
    xi = x(i);

    % k is index into xp of where we start 3 points used to form parabola

    while (xp(k) < xi) && (k < length(xp)) 
      k++;
    end
    
    %printf("xi: %f k = %d\n", xi, k);
    if k == 1
      % linear interpolation at low end
      x1 = 0; y1 = 0;
      x2 = xp(k); y2 = yp(k); 
      b = (y2-y1)/(x2-x1);
      y(i) = b*(xi-x2) + y2;
      %printf("lin1 k: %d i: %d xi: %f x1: %f y1: %f\n", k, i, xi, x1, y1);
    elseif k < length(xp)
      % parabolic interpolation
      x1 = xp(k-1); y1 = yp(k-1); 
      x2 = xp(k); y2 = yp(k); 
      x3 = xp(k+1); y3 = yp(k+1);
      a = ((y3-y2)/(x3-x2)-(y2-y1)/(x2-x1))/(x3-x1);
      b = ((y3-y2)/(x3-x2)*(x2-x1)+(y2-y1)/(x2-x1)*(x3-x2))/(x3-x1);
      y(i) = a*(xi-x2)^2 + b*(xi-x2) + y2;
      %printf("para1 k: %d i: %d xi: %f x1: %f y1: %f\n", k, i, xi, x1, y1);
    elseif (k == length(xp)) && (xi < xp(k))
      % parabolic interpolation, but shift xp points back by 1
      x1 = xp(k-2); y1 = yp(k-2); 
      x2 = xp(k-1); y2 = yp(k-1); 
      x3 = xp(k); y3 = yp(k);
      a = ((y3-y2)/(x3-x2)-(y2-y1)/(x2-x1))/(x3-x1);
      b = ((y3-y2)/(x3-x2)*(x2-x1)+(y2-y1)/(x2-x1)*(x3-x2))/(x3-x1);
      y(i) = a*(xi-x2)^2 + b*(xi-x2) + y2;
      %printf("para2 k: %d i: %d xi: %f x1: %f y1: %f\n", k, i, xi, x1, y1);
    elseif k == length(xp)
      % linear interpolation at high end
      x1 = xp(k); y1 = yp(k); 
      x2 = xp_max; y2 = 0;
      b = (y2-y1)/(x2-x1);
      y(i) = b*(xi-x1) + y1;
      %printf("lin2 k: %d i: %d xi: %f x1: %f y1: %f\n", k, i, xi, x1, y1);
    end

  end
endfunction


#{
% choose largest sample in band, idea is we care more about finding
% peaks, can handle some error in frequency. x are non linear
% (arbitrary) sampling points in kHz

function y = interp_largest(f0_Hz, AmdB, x_kHz)
  L = length(AmdB);
  x = x_kHz*1000;
  y = zeros(1,length(x));
  bw = x(2) - x(1);
  k = 1;

  for i=1:length(x)

    % determine limits of this band

    if i>1
      bw = x(i) - x(i-1);
    end
    band_low = x(i) - bw/2; band_high = x(i) + bw/2;

    % map band limits to harmonics    
 
    if x(i) < f0_Hz
      m_low = m_high = 1;
    else
      m_low = round(band_low/f0_Hz); m_high = round(band_high/f0_Hz)-1;
      m_low = max(1, m_low); m_high = min(L, m_high); m_high = max(m_low, m_high);
    end

    printf("L: %d f0: %f i: %d band_low: %f band_high: %f m_low: %d m_high: %d\n",L, f0_Hz, i, band_low, band_high, m_low, m_high);
    % find max in band

    y(i) = max(AmdB(m_low:m_high));
  end

endfunction
#}

% simple linear interpolator

function y = interp_linear(xp, yp, x)
  assert( (length(xp) == 2) && (length(yp) == 2) );

  m = (yp(2) - yp(1))/(xp(2) - xp(1));
  c = yp(1) - m*xp(1);

  y = zeros(1,length(x));
  for i=1:length(x)
    y(i) = m*x(i) + c;
  end
endfunction

% quantise input sample to nearest value in table, optionally return binary code
% should work with vectors and weights too

function [quant_out best_i bits] = quantise(levels, quant_in, weights=1)

  % if a scaler quantiser make it a col vector
  if rows(levels) == 1
    levels = levels';
  end
 
  % find closest quantiser level

  best_se = 1E32;
  for i=1:length(levels)
    se = sum( ((levels(i,:) - quant_in).^2) .* weights );
    if se < best_se
      quant_out = levels(i,:);
      best_se = se;
      best_i = i;
    end
  end

  % convert index to binary bits

  bits = index_to_bits(best_i-1, ceil(log2(length(levels))));
  
endfunction


% Quantisation functions for Wo in log freq domain

function index = encode_log_Wo(Wo, bits)
    Wo_levels = 2.^bits;
    Wo_min = 2*pi/160;
    Wo_max = 2*pi/20;

    norm = (log10(Wo) - log10(Wo_min))/(log10(Wo_max) - log10(Wo_min));
    index = floor(Wo_levels * norm + 0.5);
    index = max(index, 0);
    index = min(index, Wo_levels-1);
endfunction

function Wo = decode_log_Wo(index, bits)
    Wo_levels = 2.^bits;
    Wo_min = 2*pi/160;
    Wo_max = 2*pi/20;

    step = (log10(Wo_max) - log10(Wo_min))/Wo_levels;
    Wo   = log10(Wo_min) + step*index;

    Wo = 10 .^ Wo;
endfunction


% convert index to binary bits

function bits = index_to_bits(value, numbits)
  levels = 2.^numbits;
  bits = zeros(1, numbits);
  for b=1:numbits
    bits(b) = bitand(value,2^(numbits-b)) != 0;
  end
end

function value = bits_to_index(bits, numbits)
  value = sum(2.^(numbits-1:-1:0) .* bits);
endfunction


% Determine a phase spectra from a magnitude spectra
% from http://www.dsprelated.com/showcode/20.php
% Haven't _quite_ figured out how this works but have to start somewhere ....
%
% TODO: we may be able to sample at a lower rate, like mWo
%       but start with something that works
function [phase Gdbfk s Aw] = determine_phase(model, f, Nfft=512, ak)
  Fs      = 8000;
  max_amp = 80;
  L       = min([model(f,2) max_amp-1]);
  Wo      = model(f,1);

  sample_freqs_kHz = (Fs/1000)*[0:Nfft/2]/Nfft;           % fft frequency grid (nonneg freqs)
  Am = model(f,3:(L+3));
  AmdB = 20*log10(Am);
  rate_L_sample_freqs_kHz = (1:L)*Wo*4/pi;
  Gdbfk = interp_lanczos(rate_L_sample_freqs_kHz, AmdB, Fs/(2*1000), sample_freqs_kHz);    
 
  % Gdbfk = resample_mask(model, f, mask_sample_freqs_kHz);

  % optional input of aks for testing

  if nargin == 4
    Aw = 1 ./ fft(ak,Nfft);
    Gdbfk = 20*log10(abs(Aw(1:Nfft/2+1)));
  end

  [phase s] = mag_to_phase(Gdbfk, Nfft);

endfunction


% Non linear sampling of frequency axis, reducing the "rate" is a
% first step before VQ

function mel = ftomel(fHz)
  mel = floor(2595*log10(1+fHz/700)+0.5);
endfunction


function rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K, fstart_hz=100, fend_hz=0.95*4000)
  mel_start = ftomel(fstart_hz); mel_end = ftomel(fend_hz); 
  step = (mel_end-mel_start)/(K-1);
  mel = mel_start:step:mel_end;
  rate_K_sample_freqs_Hz = 700*((10 .^ (mel/2595)) - 1);
  rate_K_sample_freqs_kHz = rate_K_sample_freqs_Hz/1000;
endfunction


function plot_mel_sample_freqs(K, f_start_hz, f_end_hz)
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K, f_start_hz, f_end_hz);
  figure(1); clf;
  plot(rate_K_sample_freqs_kHz,'+');
endfunction


function [rate_K_surface rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model, K, Fs=8000, interp_alg='lanc') 
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K, 100, 0.95*Fs/2);
  rate_K_surface = resample_const_rate_f(model, rate_K_sample_freqs_kHz, Fs, interp_alg);
endfunction


% Resample Am from time-varying rate L=floor(pi/Wo) to fixed rate K.  This can be viewed
% as a 3D surface with time, freq, and ampitude axis.

function [rate_K_surface rate_K_sample_freqs_kHz] = resample_const_rate_f(model, rate_K_sample_freqs_kHz, Fs, interp_alg='lanc')

  % convert rate L=pi/Wo amplitude samples to fixed rate K

  max_amp = 160;
  [frames col] = size(model);
  K = length(rate_K_sample_freqs_kHz);
  rate_K_surface = zeros(frames, K);

  for f=1:frames
    Wo = model(f,1);
    L = min([model(f,2) max_amp-1]);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);

    clip_en = 0;
    if clip_en
      % clip between peak and peak -50dB, to reduce dynamic range

      AmdB_peak = max(AmdB);
      AmdB(find(AmdB < (AmdB_peak-50))) = AmdB_peak-50;
    end

    rate_L_sample_freqs_kHz = (1:L)*Wo*Fs/(2000*pi);
    
    %rate_K_surface(f,:) = interp1([0 rate_L_sample_freqs_kHz (Fs/2000)], [AmdB(1) AmdB AmdB(L)], rate_K_sample_freqs_kHz, "spline");
    if strcmp(interp_alg, 'para')
      rate_K_surface(f,:)  = interp_para(rate_L_sample_freqs_kHz, AmdB, Fs/(2*1000), rate_K_sample_freqs_kHz);
    end
    if strcmp(interp_alg, 'lanc')
      rate_K_surface(f,:)  = interp_lanczos(rate_L_sample_freqs_kHz, AmdB, Fs/(2*1000), rate_K_sample_freqs_kHz);
    end

    % equalise energy betweerate L and rate K samples.  This was required for interframe decimation/interpolation
    % when Wo is jumping about, e.g. UV frames

    energy_L = sum(Am.^2);
    rate_K_vec_lin = 10 .^ (rate_K_surface(f,:)/20);
    energy_K = sum(rate_K_vec_lin .^ 2);
    g = sqrt(energy_L/energy_K);
    rate_K_surface(f,:) += 20*log10(g);
    %printf("%d\n", f);
  end
  %printf("\n");
endfunction


% Take a rate K surface and convert back to time varying rate L

function [model_ AmdB_] = resample_rate_L(model, rate_K_surface, rate_K_sample_freqs_kHz, Fs=8000, interp_alg='lanc')
  max_amp = 160; K = columns(rate_K_surface);
  [frames col] = size(model);

  AmdB_ = zeros(frames, max_amp);
  model_ = zeros(frames, max_amp+2);
  for f=1:frames
    Wo = model(f,1);
    L = model(f,2);
    rate_L_sample_freqs_kHz = (1:L)*Wo*Fs/(2000*pi);
    
    % back down to rate L

    % dealing with end effects is an ongoing issue.....need a better solution
   
    if strcmp(interp_alg, 'para')
      AmdB_(f,1:L) = interp_para(rate_K_sample_freqs_kHz, rate_K_surface(f,:), Fs/(2*1000), rate_L_sample_freqs_kHz);    
    end
    if strcmp(interp_alg, 'lanc')
      AmdB_(f,1:L) = interp_lanczos(rate_K_sample_freqs_kHz, rate_K_surface(f,:), Fs/(2*1000), rate_L_sample_freqs_kHz);    
    end
    if strcmp(interp_alg, 'lancmel')
      rate_K_sample_freqs_mel = ftomel(rate_K_sample_freqs_kHz*1000);
      rate_L_sample_freqs_mel = ftomel(rate_L_sample_freqs_kHz*1000);
      AmdB_(f,1:L) = interp_lanczos(rate_K_sample_freqs_mel, rate_K_surface(f,:), Fs/(2*1000), rate_L_sample_freqs_mel);
    end

    % equalise energy
    
    Am_ = 10 .^ (AmdB_(f, 1:L)/20);
    energy_L = sum(Am_.^2);
    rate_K_vec_lin = 10 .^ (rate_K_surface(f,:)/20);
    energy_K = sum(rate_K_vec_lin .^ 2);
    g = sqrt(energy_K/energy_L);
    Am_ *= g;
    
    model_(f,1) = Wo; model_(f,2) = L; model_(f,3:(L+2)) = Am_;
   end
endfunction


% Decoder side interpolation of Wo and voicing, to go from 25 Hz
% sample rate used over channel to 100Hz internal sample rate of Codec
% 2.

function [Wo_ voicing_] = interp_Wo_v(Wo1, Wo2, voicing1, voicing2)
    M = 4;
    max_amp = 80;

    Wo_ = zeros(1,M); 
    voicing_ = zeros(1,M);
    if !voicing1 && !voicing2
       Wo_(1:M) = 2*pi/100;
    end

    if voicing1 && !voicing2
       Wo_(1:M/2) = Wo1;
       Wo_(M/2+1:M) = 2*pi/100;
       voicing_(1:M/2) = 1;
    end

    if !voicing1 && voicing2
       Wo_(1:M/2) = 2*pi/100;
       Wo_(M/2+1:M) = Wo2;
       voicing_(M/2+1:M) = 1;
    end

    if voicing1 && voicing2
      Wo_samples = [Wo1 Wo2];
      Wo_(1:M) = interp_linear([1 M+1], Wo_samples, 1:M);
      voicing_(1:M) = 1;
    end

    #{
    printf("f: %d f+M/2: %d Wo: %f %f (%f %%) v: %d %d \n", f, f+M/2, model(f,1), model(f+M/2,1), 100*abs(model(f,1) - model(f+M/2,1))/model(f,1), voicing(f), voicing(f+M/2));
    for i=f:f+M/2-1
      printf("  f: %d v: %d v_: %d Wo: %f Wo_: %f\n", i, voicing(i), voicing_(i), model(i,1),  model_(i,1));
    end
    #}
endfunction


function [rate_K_vec_corrected orig_error error nasty_error_log nasty_error_m_log] = correct_rate_K_vec(rate_K_vec, rate_K_sample_freqs_kHz, AmdB, AmdB_, K, Wo, L, Fs)

    % aliasing correction --------------------------------------

    % The mel sample rate decreases as frequency increases. Look for
    % any regions above 1000Hz where we have missed definition of a
    % spectral peak (formant) due to aliasing.  Adjust the rate K
    % sample levels to restore peaks.  Theory is that correct
    % definition of a formant is less important than the frequency of
    % the formant.  As long as we define a formant in that general
    % frequency area it will sound OK.

    Am_freqs_kHz = (1:L)*Wo*Fs/(2000*pi);

    % Lets see where we have made an error

    error = orig_error = AmdB(1:L) - AmdB_(1:L);

    Ncorrections = 3;      % maximum number of rate K samples to correct
    error_thresh = 3;      % only worry about errors larger than thresh

    start_m = floor(L*1000/(Fs/2));
    error(1:start_m) = 0;  % first 1000Hz is densly sampled so ignore
    nasty_error_m_log = []; nasty_error_log = [];


    rate_K_vec_corrected = rate_K_vec;
    for i=1:Ncorrections
      [mx mx_m] = max(error);

      if mx > error_thresh
        nasty_error_log = [nasty_error_log mx];
        nasty_error_m_log = [nasty_error_m_log mx_m];

        % find closest rate K sample to nasty error

        nasty_error_freq = mx_m*Wo*Fs/(2*pi*1000);
        [tmp closest_k] = min(abs(rate_K_sample_freqs_kHz - nasty_error_freq));
        rate_K_vec_corrected(closest_k) = AmdB(mx_m);

        % zero out error in this region and look for another large error region

        k = max(1, closest_k-1); 
        rate_K_prev_sample_kHz = rate_K_sample_freqs_kHz(k);
        k = min(K, closest_k+1); 
        rate_K_next_sample_kHz = rate_K_sample_freqs_kHz(k);

        [tmp st_m] = min(abs(Am_freqs_kHz - rate_K_prev_sample_kHz));
        [tmp en_m] = min(abs(Am_freqs_kHz - rate_K_next_sample_kHz));
        if closest_k == K
         en_m = L;
        end 
        error(st_m:en_m) = 0;
      end
    end
endfunction


% Given a vector of rate K samples, deltaf encodes/decodes delta
% amplitude using a huffman encoder, returning quantised samples and
% bit stream

function [rate_K_vec_ bits] = deltaf_quantise_rate_K_huff(rate_K_vec, E, nbits_max)
  K = length(rate_K_vec);
  nbits_remaining = nbits_max;
  
  % start with k=3, around 250Hz, we assume that's quantised as the
  % mean frame energy, as samples before that might be stuck in the
  % HPF.

  rate_K_vec_ = zeros(1,K);
  rate_K_vec_(3) = E;
  
  % encoding of differences
  
  levels = [0 6 -6 -12 12]; symbols = {[0 0],[1 0],[1 1],[0 1 0;],[0 1 1]};

  % this is pretty coarse (note no 0dB level) but sounds OK
  
  %levels = [-6 +6 -12 +12]; symbols = {[0 0],[1 0],[1 1],[0 1]};

  % move backwards to get target for first two samples
  
  bits = [];
  for m=2:-1:1
    target = rate_K_vec(m) - rate_K_vec_(m+1);
    [target_ best_i] = quantise(levels, target);
    bits = [bits symbols{best_i}];
    nbits_remaining -= length(symbols{best_i});
    rate_K_vec_(m) = rate_K_vec_(m+1) + target_;
  end
  
  % then forwards for rest of the target samples
  
  for m=4:K
    target = rate_K_vec(m) - rate_K_vec_(m-1);
    if nbits_remaining >= 3 
      [target_ best_i] = quantise(levels, target);
      bits = [bits symbols{best_i}];
      nbits_remaining -= length(symbols{best_i});
    elseif nbits_remaining == 2
      [target_ best_i] = quantise(levels(1:3), target);
      bits = [bits symbols{best_i}];
      nbits_remaining = 0;
    elseif nbits_remaining < 2
      target_ = -6; % if we've run out of bits just tail off
    end
    rate_K_vec_(m) = target_ + rate_K_vec_(m-1);    
    %printf("m: %d length: %d nbits_remaining: %d\n", m, length(bits), nbits_remaining);
  end

  bits = [bits zeros(1,nbits_remaining)];
endfunction


% Given a vector of bits, and the k=3 frame amplitude sample E, decode to a
% rate K vector of huffman encoded quantised spectral amplitude samples

function rate_K_vec_ = deltaf_decode_rate_K_huff(bits, E, K)
  
  % start with k=3, around 250Hz, we assume that's quantised as the
  % mean frame energy, as samples before that might be stuck in the
  % HPF.

  rate_K_vec_ = zeros(1,K);
  rate_K_vec_(3) = E;
  
  % move backwards to get target for first two samples
  
  for m=2:-1:1
    [target_ bits] = deltaf_dec_one_symbol(bits);
    rate_K_vec_(m) = rate_K_vec_(m+1) + target_;
  end
  
  % then forwards for rest of the target samples
  
  for m=4:K
    if length(bits) >= 2
      [target_ bits] = deltaf_dec_one_symbol(bits);
    else
      target_ = -6; % if we've run out of bits just tail off
    end
    rate_K_vec_(m) = target_ + rate_K_vec_(m-1);    
    %printf("m: %d nbits_remaining: %d\n", m, length(bits));
  end
endfunction


% decode one symbol and truncate bits array

function [target_ bits] = deltaf_dec_one_symbol(bits)
  levels = [0 6 -6 -12 12]; symbols = {[0 0],[1 0],[1 1],[0 1 0;],[0 1 1]};

  two_bits = bits(1:2);

  if two_bits == [0 1]
    if length(bits) > 2
      % OK a three bit code
      if bits(3) == 0
        target_ = levels(4);
      else
        target_ = levels(5);
      end
      bits = bits(4:end);
    else
      # we must have a bit error
      target_ = levels(5);
      bits = [];
    end  
  else
    if two_bits == [0 0]
      target_ = levels(1);
    end
    if two_bits == [1 0]
      target_ = levels(2);
    end
    if two_bits == [1 1]
      target_ = levels(3);
    end      
    bits = bits(3:end);
  end
endfunction


function rate_K_vec_ = dct_quantise_rate_K(rate_K_vec)
    K = length(rate_K_vec);
    D = dct(rate_K_vec);
    printf("\n");
    D_ = 8*round(D/8);
    D_(1) = D(1);
    for d=1:K
      printf("%4d",round(D(d)/8));
    end
    rate_K_vec_ = idct(D_);
endfunction


function un(sl)
  U = unique(sl,"rows")
  cnt = [];
  for v=1:length(U)
    sm = sum(ismember(sl,U(v,:),"rows"));
    cnt = [cnt sm];
  end
  s = sort(cnt, "descend");
  figure(1); clf; subplot(211); plot(s); subplot(212); plot(cumsum(s));  
endfunction


% Joint Wo and LPC energy vector quantiser developed by Jean-Marc Valin.
% Octave port of functions in quantise.c

function w = compute_weights2(x, xp)
  w(1) = 30;
  w(2) = 1;
  if x(2) < 0
     w(1) *= 0.6;
     w(2) *= 0.3;
  end
  if x(2) < -10
     w(1) *= 0.3;
     w(2) *= 0.3;
  end

  % Higher weight if pitch is stable

  if abs(x(1)-xp(1)) < 0.2
     w(1) *= 2;
     w(2) *= 1.5;
  elseif abs(x(1)-xp(1)) > 0.5
     % Lower if not stable
     w(1) *= 0.5;
  end
  
  % Lower weight for low energy
  
  if x(2) < xp(2) - 10
     w(2) *= 0.5;
  end
  if x(2) < xp(2) - 20
     w(2) *= .5;
  end

  % Square the weights because it's applied on the squared error
  
  w(1) *= w(1);
  w(2) *= w(2);
endfunction


function [Wo_ E_ xq] = quantise_WoE(Wo, E, xq, vq)
  ge_coeff = [0.8 0.9];

  % VQ is only trained for Fs = 8000 Hz

  Fs = 8000;            
  Fo_min = 50; Fo_max = 400;
  P_min = Fs/Fo_max; P_max = Fs/Fo_min;
  Wo_min = 2*pi/P_max;
  Wo_max = 2*pi/P_min;

  E = max(1,E);
  
  x(1) = log10(Wo/Wo_min)/log10(2);
  x(2) = 10.0*log10(1e-4 + E);

  w = compute_weights2(x, xq);
  w = [30 1];
  err  = x - ge_coeff .* xq;
  err_ = quantise(vq, err, w);
  %err_ =  err;
  xq = ge_coeff .* xq + err_;

  #{
    x = log2(4000*Wo/(PI*50));
    2^x = 4000*Wo/(PI*50)
    Wo = (2^x)*(PI*50)/4000;
  #}

  Wo_ = (2 ^ xq(1))*Wo_min;

  Wo_ = min(Wo_max,Wo_);
  Wo_ = max(Wo_min,Wo_);

  E_ = 10.0 ^ (xq(2)/10.0);
  E_ = max(1,E_);
endfunction


% Find distance of each vector in codebook using gain shape search
% Uses discrete values, in unit range

function dist_gain_shape(codebook, target)
  [entries K] = size(codebook);
  dist = zeros(entries,1);
  best = zeros(entries,K);
  for i=1:entries
    min_dist = 1000;
    for g=-10:20
       adist = sum(abs(target - (codebook(i,:) + g)));
       if adist < min_dist
         min_dist = adist;
         min_g    = g;
       end
    end
    dist(i) = min_dist;
    best(i,:) = codebook(i,:) + min_g;
  end
  figure(1); clf; plot(dist);
  axis([1 entries 0 K])
  target
  [tmp min_i] = min(dist);
  target - best(min_i,:)
endfunction


% As above but without the gain

function dist_shape(codebook, target)
  [entries K] = size(codebook);
  dist = zeros(entries,1);
  for i=1:entries
    dist(i) = sum(abs(target - (codebook(i,:))));
  end
  figure(1); clf; plot(dist);
  axis([1 entries 0 K])
  target
  [tmp min_i] = min(dist);
  target - codebook(min_i,:)
endfunction


% Given a vector of rate K samples, deltaf encodes/decodes delta
% amplitude using a fixed bit/sample allocation, returning quantised
% samples and bit stream

function [rate_K_vec_ bits] = deltaf_quantise_rate_K_fixed(rate_K_vec, E, nbits_max)
  K = length(rate_K_vec);
  
  % start with k=3, around 250Hz, we assume that's quantised as the
  % mean frame energy, as samples before that might be stuck in the
  % HPF.

  rate_K_vec_ = zeros(1,K);
  rate_K_vec_(3) = E;
  
  % encoding of differences
  
  levels_2bit =  [-3 +3 -9 9];
  levels_3bit =  [0 -6 +6 -12 +12 -18 +18 -24];

  % move backwards to get target for first two samples
  
  bits = [];
  for m=2:-1:1
    target = rate_K_vec(m) - rate_K_vec_(m+1);
    [target_ best_i] = quantise(levels_2bit, target);
    % printf("m: %d target_ %f best_i: %d\n", m, target_, best_i);
    bits = [bits index_to_bits(best_i-1, 2)];
    rate_K_vec_(m) = rate_K_vec_(m+1) + target_;
  end
  
  % then forwards for rest of the target samples
  
  for m=4:9
    target = rate_K_vec(m) - rate_K_vec_(m-1);
    [target_ best_i] = quantise(levels_3bit, target);
    bits = [bits index_to_bits(best_i-1, 3)];
    rate_K_vec_(m) = target_ + rate_K_vec_(m-1);    
  end
  for m=10:K
    target = rate_K_vec(m) - rate_K_vec_(m-1);
    [target_ best_i] = quantise(levels_2bit, target);
    bits = [bits index_to_bits(best_i-1, 2)];
    rate_K_vec_(m) = target_ + rate_K_vec_(m-1);    
  end
  assert(length(bits) == nbits_max);
endfunction


% Given a vector of bits, and the k=3 frame amplitude sample E, decode
% to a rate K vector of encoded quantised spectral amplitude samples
% with fixed bit allocation

function rate_K_vec_ = deltaf_decode_rate_K_fixed(bits, E, K)
  
  % start with k=3, around 250Hz, we assume that's quantised as the
  % mean frame energy, as samples before that might be stuck in the
  % HPF.

  rate_K_vec_ = zeros(1,K);
  rate_K_vec_(3) = E;
  
  levels_2bit = [-3 +3 -9 9];
  levels_3bit =  [0 -6 +6 -12 +12 -18 +18 -24];

  % move backwards to get target for first two samples
  
  for m=2:-1:1
    index = bits_to_index(bits(1:2),2) + 1;
    target_ = levels_2bit(index);
    rate_K_vec_(m) = rate_K_vec_(m+1) + target_;
    bits = bits(3:end);
  end
  
  % then forwards for rest of the target samples
  
  for m=4:9
    index = bits_to_index(bits(1:3), 3) + 1;
    target_ = levels_3bit(index);
    rate_K_vec_(m) = rate_K_vec_(m-1) + target_;
    bits = bits(4:end);
  end
  for m=10:K
    index = bits_to_index(bits(1:2), 2) + 1;
    target_ = levels_2bit(index);
    rate_K_vec_(m) = rate_K_vec_(m-1) + target_;
    if m != K
      bits = bits(3:end);
    end
  end
  
endfunction


#{
  DCT coeff scalar quantiser, analysis stage.  Plots PDFs, returns a
  vector of quantiser levels.

  Input is matrix of K columns row-vectors that represent one frame of
  spectrum at rate K.
#}

function quantiser_levels = build_dct_quantiser(train_surf, qstepdB=1)
  [nr K] = size(train_surf);

  % remove low energy rows

  m = mean(train_surf');
  figure(4); plot(m);
  ind = find(m>10);
  train_surf = train_surf(ind,:);
  nr2 = length(ind);
  %nr2 = nr;
  printf("K: %d nr: %d nr2: %d\n", K, nr, nr2);

  D = dct(train_surf')';

  figure(1); clf;
  plot(std(D));
  title('Std Dev of each DCT coeff');

  figure(2); clf;
  nr = ceil(sqrt(K)); nc = ceil(K/nr);

  Tbits = 0; quantiser_levels = [];
   
  for k=1:K
    subplot(nr, nc, k)
    v = D(:,k);
    printf("k: %d mean %5.2f std: %5.2f min: %5.2f max: %5.2f\n", k, mean(v), std(v), min(v), max(v));
    q_max = mean(v)+2*std(v); q_min = mean(v)-2*std(v);
    levels = q_min:qstepdB:q_max;
    nlevels = length(levels);
    Tbits += log2(nlevels);
    printf("    quantiser: min: %4.2f max: %4.2f nlevels: %d bits: %2.1f\n", q_min, q_max, nlevels, log2(nlevels));

    % limit for histogram
    
    v = min(v, q_max);
    v = max(v, q_min);
    [nn xx] = hist(v,50);
    bar (xx, nn)

    v_ = zeros(nr2,1);
    for r=1:nr2
      v_(r) = quantise(levels, v(r));
    end
    
    E(:,k) = v_;

    quantiser_levels = [quantiser_levels; q_min q_max]; 
  end

  train_surf_ = idct(E')';
  error = train_surf_ - train_surf;
  mse = mean(mean(error .^ 2));
  figure(3)
  mesh(error(1:1000,:))
  printf("mse: %4.2f Tbits: %d\n", mse, Tbits);
endfunction


function bits = bits_for_this_symbol(symbols, s)
  % bits/sym for top 5, assume rest 5 bits
  %      11 10 00 010 0110 0111
  bps = [ 2  2  2   3    4    4]; 
  
  ind = find(symbols == s);
  if ind <= length(bps)
    bits = bps(ind);
  else
    bits = 5;
  end
endfunction


function [surf_ D E] = dct_quantiser(surf, quantiser_levels, method=1, qstepdB=1, limit=0)
  [nr K] = size(surf);

  printf("K: %d nr: %d\n", K, nr);

  % clamp lower limit of mean to 10dB
  #{
  m = mean(surf')';
  surf -= m;
  m = min(10,m);
  surf += m;
  #}
  
  D = dct(surf')';

  Tbits = 0;
   
  if method == 1
    for k=1:K
      q_min = quantiser_levels(k,1); q_max = quantiser_levels(k,2);
      levels = q_min:qstepdB:q_max;
      nlevels = length(levels);
      if nlevels == 1
        levels = (q_max + q_min)/2;
      end
      if nlevels == 2
        m = (q_max + q_min)/2;
        levels = [(m - qstepdB/2) (m + qstepdB/2)];
      end
      if nlevels == 3
        m = (q_max + q_min)/2;
        levels = [q_min m q_max];
      end
      Tbits += log2(nlevels);
      printf("    quantiser: min: %4.2f max: %4.2f nlevels: %d bits: %2.1f\n", q_min, q_max, nlevels, log2(nlevels));

      v_ = zeros(nr,1);
      v = D(:,k);
      for r=1:nr
        v_(r) = quantise(levels, v(r));
      end
    
      E(:,k) = v_;
    end
  end
  
  if method == 2
    % quantise
    E = round(D/qstepdB);
    % count symbols
    symbols = []; count = [];
    [nr nc]= size(E);
    if limit
      nc = limit-1;
    end
    for r=1:nr
      for c=2:nc
        s = E(r,c);
        ind = find(symbols == s);
        if length(ind)
          count(ind)++;
        else
          symbols = [symbols s];
          count(length(symbols)) = 1;
        end
      end
    end

    % sort into order

    [count ind] = sort(count, "descend");
    symbols = symbols(ind);

    % estimate bits/symbol by huffman coding direct and differences.  Clever part is we
    % choose method based on min bits, which I guess costs an extra bit.  It's clever
    % because we need direct for quick transitions, but during steady voiced parts the
    % changes are small so coding differences does a good job.

    Tbits = 0; Nsyms = 0;
    Nbits_direct_log = Nbits_diff_log = Nbits_log = zeros(1,nr);
   
    for r=3:nr
      Nbits_direct = 0; Nbits_diff = 0;
      for c=2:nc
        s_direct = E(r,c);
        s_diff = E(r,c) - E(r-2,c);
        Nbits_direct += bits_for_this_symbol(symbols, s_direct);
        Nbits_diff   += bits_for_this_symbol(symbols, s_diff);
        Nsyms++;
      end
      Nbits_direct_log(r) += Nbits_direct;
      Nbits_diff_log(r) += Nbits_diff;
      Nbits_log(r) += min(Nbits_direct, Nbits_diff);

      Tbits += min(Nbits_direct, Nbits_diff);
    end
    figure(5); clf; plot(Nbits_direct_log(3:end),'b');
    smooth4 = conv(Nbits_log(3:end),[1 0 1 0 1 0 1 0])/4;
    hold on; plot(Nbits_diff_log(3:end),'g'); plot(Nbits_log(3:end),'r'); hold off;
    figure(6); clf; plot(smooth4,'m');
    
    ent = 0;
    for i=1:length(symbols)
      wi = count(i)/sum(count);
      printf("%2d %4d %6d %4.3f %4.3f\n", i, symbols(i), count(i), wi, -wi*log2(wi));
      ent += -wi*log2(wi);
    end
    printf("mean bits/frame: %3.1f mean bits/sym: %3.2f entropy: %3.2f bits\n",  mean(Nbits_log), Tbits/Nsyms, ent);
    if limit
      E(:, limit:end) = 0;
    end
    E *= qstepdB;
  end

  surf_ = idct(E')';
  error = surf_ - surf;
  mse = mean(mean(error .^ 2));
  figure(3)
  [nr nc] = size(error);
  nr = min(nr,100);
  mesh(error(1:nr,:))
  figure(4)
  mesh(D(25:50,1:10))
  printf("mse: %4.2f Tbits: %d\n", mse, Tbits);
endfunction


# Generate a huffman code from a matrix (surface) of spectral
# magnitudes.  The Huffman code can be used for encoding the DCT of the
# rows of the surface (mag samples of each frame).
#
#   Set qstepdB to the quanisation step size, e.g. 3 or 6dB is about where we can
#   notice the effect of quantisation of the DCTs (and spectral magnitides)
#
#   Set "max_dcts" to max number of dcts coeffs you will quantise, as this affects
#   probabilities (we get many zeros in high order coeffs)
#
#   octave:49> newamp2; p_table = design_huffman_enc(all_surf(:,2:35), 6, 18);

function [symbols huff] = design_huffman_enc(surf, qstepdB=6, max_dcts=18)
  [nr K] = size(surf);

  printf("K: %d nr: %d qstepdB: %3.2f\n", K, nr, qstepdB);

  D = dct(surf')';

  % quantise to step size in dB
    
  E = round(D/qstepdB);

  % count symbols, ignoring first (DC) coeff which we will scalar quantise

  symbols = []; count = [];
  [nr nc]= size(E);
  if max_dcts
    nc = max_dcts+1;
    printf("cols 2 to %d (%d total)\n", nc, nc-2+1);
  end
  for r=1:nr
    for c=2:nc
      s = E(r,c);
      ind = find(symbols == s);
      if length(ind)
        count(ind)++;
      else
        symbols = [symbols s];
        count(length(symbols)) = 1;
      end
    end
  end

  % sort into order

  [count ind] = sort(count, "descend");
  symbols = symbols(ind);

  Nsymbols = sum(count);
  printf("Nsymbols = %d\n", Nsymbols);
  
  % estimate entropy
    
  H = 0;
  p_table = [];
  printf(" i symb  count  prob    wi\n");
  for i=1:length(symbols)
    wi = count(i)/Nsymbols; p_table = [p_table wi];
    printf("%2d %4d %6d %4.3f %4.3f\n", i, symbols(i), count(i), wi, -wi*log2(wi));
    H += -wi*log2(wi);
  end

  % design Huffman code

  huff = huffmandict (1, p_table, 1);
  L = 0;
  for i=1:length(huff)
    L += p_table(i)*length(huff{i});
  end

  printf("Entropy: %3.2f bits/symbol  Huffman code: %3.2f bits/symbol\n",  H, L);
endfunction


# Huffman encodes (and decodes) a symbol, if input symbols is out of
# range of quantiser we choose nearest symbol

function [s_ bits] = huffman_enc_symb(symbols, huff, s)
  min_dist = 1E32; ind = 1;
  for i=1:length(symbols)
    dist = (symbols(i) - s) .^ 2;
    if dist < min_dist
      ind = i;
      min_dist = dist;
    end
  end

  s_ = symbols(ind);
  bits = huff{ind};
endfunction


% Huffman decode a bit stream of symbols.  Terminates list if we get a
% bit error in decode

function [s_ error_flag] = huffman_decode_bitstream(symbols, huff, bits)
  error_flag = 0;
  min_cw_length = length(huff{1});
  match = 1;
  s_ = [];
  
  while ((length(bits) >= min_cw_length) && match)

    % search through list of codes to find a match

    match = 0;
    for i=1:length(symbols)
      cw = huff{i}; lcw = length(cw);
      if (length(bits) >= lcw) & !match
        match = isequal(cw, bits(1:lcw));
        if match
          s_ = [s_ symbols(i)];
          bits = bits(lcw+1:end);
        end
      end
    end

    % if no match found, say due to bit error, we drop out of loop
  end

  if (length(bits) > min_cw_length) && !match
    error_flag = 1;
  end
endfunction

 
% Huffman encodes (and decodes) the DCTS of a surface, except first (DCT) coeff

function [surf_ dc bits_surf] = huffman_encode_surf(surf, qstepdB=6, max_dcts=18, max_bits=100, symbols, huff)
  dec = 2;
  [nr K] = size(surf);

  % allow room for direct/diff bit
  
  max_bits_huff = max_bits - 1;
  
  printf("K: %d nr: %d qstepdB: %3.2f max_dcts: %d\n", K, nr, qstepdB, max_dcts);

  % limit num DCTs we encode (nc) to less than K to save bits, high
  % order DCTs tend to be small

  nc = K;
  if max_dcts
    nc = max_dcts+1;
    printf("cols 2 to %d (%d total)\n", nc, nc-2+1);
  end

  % DCT and initial quantisation to step size
  
  D = dct(surf')';
  E = D/qstepdB;
  dc = E(:,1);
  
  % bit stream for each row (frame) is stored in a cell array
  
  bits_surf = cell(nr,1);
  Nbits_log = Nbits_direct_log = Nbits_diff_log = zeros(1,nr);
  E_ = zeros(nr,K); prev_row = zeros(1,nc);
  Tbits = Nsyms = 0;
  E_dec = zeros(nr,K);
  
  % encode each row

  for r=1:dec:nr
    bits_direct_row = [];  bits_diff_row = []; Nbits_direct = 0; Nbits_diff = 0;

    % DC just copied directly, quantised externally
    
    E_(r,1) = E(r,1); E_dec(r,1) = E(r,1);
    
    direct_row_ = diff_row_ = zeros(1,nc);
    ndir_row = ndiff_row = 0; len_bits_direct_row = len_bits_diff_row = 0;
    
    for c=2:nc

      % try direct quantisation

      if len_bits_direct_row < max_bits_huff
        s_direct = E(r,c);
        [s_direct_ bits_direct] = huffman_enc_symb(symbols, huff, s_direct);
        if len_bits_direct_row + length(bits_direct) <= max_bits_huff
          % can we squeeze in bits for latest symbol?
          bits_direct_row = [bits_direct_row bits_direct];
          direct_row_(c) = s_direct_;
          ndir_row++;
          len_bits_direct_row = len_bits_direct_row + length(bits_direct);
        else
          % can't fit any more symbols? Then signal we are finished
          len_bits_direct_row = max_bits_huff;
        end
      end

      % try differential quantisation
      
      if len_bits_diff_row < max_bits
        s_diff = E(r,c) - prev_row(c);
        [s_diff_ bits_diff] = huffman_enc_symb(symbols, huff, s_diff);
        if len_bits_diff_row + length(bits_diff) <= max_bits_huff
          bits_diff_row = [bits_diff_row bits_diff];
          diff_row_(c) = s_diff_;
          ndiff_row++;
          len_bits_diff_row = len_bits_diff_row + length(bits_diff);
        else
          len_bits_diff_row = max_bits_huff;
        end
      end
      
      Nsyms++;
    end

    % choose quant method with the least number of bits

    Nbits_direct = length(bits_direct_row); Nbits_diff = length(bits_diff_row);
    Nbits_direct_log(r) = Nbits_direct; Nbits_diff_log(r) = Nbits_diff;
    e_direct = sum((direct_row_(2:nc) - E(r,2:nc)) .^ 2);
    e_diff = sum((diff_row_(2:nc) + prev_row(2:nc) - E(r,2:nc)) .^ 2);
    
    %if Nbits_direct < Nbits_diff
    if e_direct <= e_diff
      bits_surf{r} = [1 bits_direct_row];
      Tbits += Nbits_direct;
      E_(r,2:nc) = direct_row_(2:nc);
      Nbits_log(r) = Nbits_direct; diff_flag(r) = 0;
    else
      bits_surf{r} = [0 bits_diff_row];
      Tbits += Nbits_diff;
      E_(r,2:nc) = diff_row_(2:nc) + prev_row(2:nc);
      Nbits_log(r) = Nbits_diff; diff_flag(r) = 1;
    end

    % pad out to max_bits
    
    bits_surf{r} = [bits_surf{r} zeros(1,max_bits - length(bits_surf{r}))];

    % test huffman bitstream decoder

    bits_dec = bits_surf{r};
    [s_dec error_flag] = huffman_decode_bitstream(symbols, huff, bits_dec(2:end));
    % printf("r: %d bits_dec(1): %d l: %d error_flag: %d\n", r, bits_dec(1), length(s_dec), error_flag);
    row_dec = zeros(1,nc);
    row_dec(2:length(s_dec)+1) = s_dec;
    if bits_dec(1)
      E_dec(r,2:nc) = row_dec(2:nc);
    else      
      E_dec(r,2:nc) = row_dec(2:nc) + prev_row(2:nc);
    end
    assert(E_dec(r,2:nc) == E_(r,2:nc));
    
    % update memory - note we use quantised symbols as that's what we have at decoder
    
    if r >=3
      prev_row = E_(r,1:nc);

      % if we are decimating, interpolate DCTs to get original frame rate
      
      if dec == 2
        E_(r-1,:) = 0.5*E_(r-2,:) + 0.5*E_(r,:);
      end
      
    end

  end

  % transform back to surface and calculate MSE

  E_ *= qstepdB;
  surf_ = idct(E_')';  
  
  error = surf_ - surf;
  mse = mean(mean(error(1:dec:end,:) .^ 2));
  
  figure(1); clf;
  [nr nc] = size(error);
  nr = min(nr,300);
  mesh(error(1:dec:nr,:))

  figure(2);
  subplot(122,"position",[0.7 0.05 0.25 0.85])
  hist(mean(error(1:dec:end,:).^2,2));
  subplot(121,"position",[0.1 0.05 0.5 0.85])
  plot(mean(error(1:dec:end,:).^2,2));
  title('Mean squared error per frame');
  
  figure(3);
  subplot(122,"position",[0.7 0.05 0.25 0.9])
  hist(Nbits_log(1:dec:end));
  subplot(121,"position",[0.1 0.05 0.5 0.9])
  plot(diff_flag(1:dec:end)*10,'b;diff flag;'); hold on; plot(Nbits_log(1:dec:end),'g;Nbits/fr;'); hold off;
  
  printf("mse: %4.2f dB^2 mean bits/frame: %3.1f mean bits/sym: %3.1f\n", mse, mean(Nbits_log(1:dec:end)), Tbits/Nsyms);
endfunction


% Returns a quantised surface from matrix of bist streams for each frame

function surf_ = huffman_decode_surf(K=34, qstepdB=6, max_dcts=18, symbols, huff, bits_surf, dc)
  dec = 2;
  [nr tmp] = size(bits_surf);
 
  printf("K: %d nr: %d qstepdB: %3.2f max_dcts: %d\n", K, nr, qstepdB, max_dcts);

  % limit num DCTs we encode (nc) to less than K to save bits, high
  % order DCTs tend to be small

  nc = K;
  if max_dcts
    nc = max_dcts+1;
    printf("cols 2 to %d (%d total)\n", nc, nc-2+1);
  end
  
  prev_row = zeros(1,nc);
  Tbits = Nsyms = 0;
  E_dec = zeros(nr,K);
  
  % decode each row

  for r=1:dec:nr
    bits_direct_row = [];  bits_diff_row = []; Nbits_direct = 0; Nbits_diff = 0;

    % DC is quantised externally
    
    E_dec(r,1) = dc(r);
    
    % huffman bitstream decoder

    bits_dec = bits_surf{r};
    [s_dec error_flag] = huffman_decode_bitstream(symbols, huff, bits_dec(2:end));
    row_dec = zeros(1,nc);
    row_dec(2:length(s_dec)+1) = s_dec;
    if bits_dec(1)
      E_dec(r,2:nc) = row_dec(2:nc);
    else      
      E_dec(r,2:nc) = row_dec(2:nc) + prev_row(2:nc);
    end
    
    % update memory for diff decoder
    
    if r >=3
      prev_row = E_dec(r,1:nc);

      % if we are decimating, interpolate DCTs to get original frame rate
      
      if dec == 2
        E_dec(r-1,:) = 0.5*E_dec(r-2,:) + 0.5*E_dec(r,:);
      end      
    end

  end

  % transform back to surface and calculate MSE

  E_dec *= qstepdB;
  surf_ = idct(E_dec')';    
endfunction

