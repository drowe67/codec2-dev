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
% Code here to support a bunch of experimental ideas, many that didn't work out.

1;
melvq; % mbest VQ functions

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

function [quant_out best_i bits] = quantise(levels, quant_in)

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
  value = 2.^(numbits-1:-1:0) * bits;
endfunction


% helper function to find polynomial coeffs for a parabola

function b = parabola_coeffs(x, y)
  A = [(x.^2)' x' [1 1 1]']
  b = inv(A)*y';    
endfunction


% generate a zig-zag linear to square mapping matrix

function map = create_zigzag_map(nr,nc)
  map = zeros(nr, nc);

  state = 'zig';
  r = c = 1;

  for i=1:nr*nc

    printf("%s r: %d c: %d i %d\n", state, r, c, i);
    map(r,c) = i;

    next_state = state;
    if state == 'zig'
      % move SE
      c -= 1; r += 1;
      if r > nr
        r = nr; c+=2;
        next_state = 'zag';
      end
      if c < 1
        c = 1;
        next_state = 'zag';
      end
    end

    if state == 'zag'
      % move SE
      r -= 1; c +=1;
      if c > nc
        c = nc; r+=2;
        next_state = 'zig';
      end
      if r < 1
        r = 1;
        next_state = 'zig';
      end
    end
    state = next_state;
  end
endfunction


% reshape matrix m as a vector v by reading out elements in zig-zag pattern

function v = mtov_zigzag(m)
  [nr nc] = size(m);
  map = zigzag(nr,nc)
  v = zeros(1,nr*nc);
  for r=1:nr
    for c=1:nc
      v(map(r,c)) = m(r,c);
    end
  end
endfunction


% extracts DCT information for rate K surface

function unwrapped_dcts = dct_blocks(surf, Nt=16)
  [frames K] = size(surf);

  % break into 160ms blocks, 2D DCT, truncate, IDCT

  Nblocks = floor(frames/Nt);
  unwrapped_dcts = zeros(Nblocks,Nt*K);

  for n=1:Nblocks
    st = (n-1)*Nt+1; en = st + Nt - 1;
    D = dct2(surf(st:en,:));
    unwrapped_dcts(n,:) = reshape(D',1,Nt*K);
  end
endfunction


% Determines a map for quantising 2D DCT coeffs in order of rms value
% (ie std dev) of each coeff.  Those coeffs with the greatest change
% need the most bits to quantise

function [map rms_map mx mx_ind unwrapped_dcts] = create_map_rms(rate_K_surface, nr, nc)
  unwrapped_dcts = dct_blocks(rate_K_surface, nr);
  [mx mx_ind] = sort(std(unwrapped_dcts));
  mx_ind = fliplr(mx_ind); mx = fliplr(mx);
  map = rms_map = zeros(nr,nc);
  for i=1:nr*nc
    r = floor((mx_ind(i)-1)/nc) + 1;
    c = mx_ind(i) - (r-1)*nc;
    %printf("%d %d %d\n", i, r, c);
    map(r,c) = i;
    rms_map(r,c) = mx(i);
  end
endfunction


% plot histogram of each 2D DCT coeff, so we can get a feel for
% quantiser design

function plot_dct2_hists(rate_K_surface, nr, nc)
  [map rms_map mx mx_ind unwrapped_dcts] = create_map_rms(rate_K_surface, nr, nc);
  Ncoeff = nr*nc;
  fign = 1; subplotn = 1;
  close all; figure(fign); clf;
  Nplot = 60;
  for i=1:Nplot   
    subplot(5,4,subplotn);
    d = unwrapped_dcts(:,mx_ind(i));
    d = round(d/4);
    hist(d);
    subplotn++;
    if (subplotn > 20) && (i != Nplot)
      subplotn = 1;
      fign++;
      figure(fign); clf;
    end
  end
endfunction


% Gather run length data for each 2D DCT coeff, to see if run length encoding
% can help

function [run_length d]= plot_run_length(rate_K_surface, nr, nc)
  [map rms_map mx mx_ind unwrapped_dcts] = create_map_rms(rate_K_surface, nr, nc);
  Ncoeff = nr*nc;
  [Nblocks tmp] = size(unwrapped_dcts);

  % first get histogram of DCT values -----------------------------------

  % some mild quantisation

  unwrapped_dcts = round(unwrapped_dcts/4);

  % note we only consider first half of DCT coeffs, unlikely to use all 

  d = [];
  for i=1:Nblocks
    d = [d unwrapped_dcts(i,mx_ind(1:Ncoeff/2))]; 
  end

  % note we remove outliers from plot as very low prob symbols

  d = d(find(abs(d)<10));
  figure(1); clf; [Wi, ii] = hist(d,-10:10,1); plot(ii,Wi);
  length(d)
  Wi = Wi(find(Wi > 0));
  %sum(Wi)
  %-log2(Wi)
  %-Wi .* log2(Wi)
  printf("bits/symbol: %2.2f\n", sum(-Wi .* log2(Wi)));

  % now measure run lengths --------------------------------------------
  
  run_length = zeros(21,Ncoeff);
  state = 'idle';

  for i=2:length(d)

    next_state = state;

    if state == 'idle'
      if d(i-1) == d(i)
        next_state = 'trac';
        arun_length = 2;
      else
        run_length(d(i)+10, 1)++;
      end
    end

    if state == 'trac'
      if d(i-1) == d(i)
        arun_length++;
      else
        next_state = 'idle';
        run_length(d(i-1)+10, arun_length)++;
      end
    end

    state = next_state;

  end

  figure(2); clf; mesh(run_length(:,1:10));
endfunction


% Design kmeans quantisers for each DCT coeff.  This didn't work very well.

function [quantisers nbits] = design_quantisters_kmeans(rate_K_surface, nr, nc, nbits_max)
  [map rms_map mx unwrapped_dcts] = create_map_rms(rate_K_surface, nr, nc);
  nq = nr*nc;
  quantisers = zeros(nq, 2^nbits_max); nbits = zeros(nq,1);
  for i=1:nq

    % work out number of levels for this quantiser such that it is a
    % power of 2 for integer number of bits

    nlevels = (2^nbits_max);
    nbits(i) = round(log2(nlevels));
    nlevels = 2 .^ nbits(i);

    if i <= 100
      printf("%d %d\n", i, nlevels);
      [idx, centers] = kmeans(unwrapped_dcts(:,i), nlevels);
      quantisers(i,1:nlevels) = sort(centers);
    end
  end
endfunction


% Uniform quantisers designed to fit limits of each DCT coeff

function [quantisers nlevels] = design_quantisters_uniform(rate_K_surface, nr, nc, nlevels_max)
  [map rms_map mx mx_ind unwrapped_dcts] = create_map_rms(rate_K_surface, nr, nc);
  
  nq = nr*nc;
  quantisers = zeros(nq, nlevels_max); nlevels = zeros(nq, 1);

  for i=1:nq
    d = unwrapped_dcts(:,mx_ind(i));
    d = floor(d/16);
    q_min = floor(min(d));
    q_max = ceil(max(d));
    nlevels(i) = q_max-q_min+1;
    quantisers(i,1:nlevels(i)) = 16*(q_min:q_max);
  end
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
  Am = model(f,3:(L+2));
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

    %printf("%d\n", f);
  end
  %printf("\n");
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
    
#{
    if pad_end
      AmdB_(f,1:L) = interp1([0 rate_K_sample_freqs_kHz Fs/2000], 
                             [rate_K_surface(f,1) rate_K_surface(f,:) rate_K_surface(f,K)], 
                             rate_L_sample_freqs_kHz, 
                             "spline");
    else
      AmdB_(f,1:L) = interp1([0 rate_K_sample_freqs_kHz], 
                             [rate_K_surface(f,1) rate_K_surface(f,:)], 
                             rate_L_sample_freqs_kHz, 
                             "spline");
    end
#}

    %AmdB_(f,1:L) = interp_para(rate_K_sample_freqs_kHz, rate_K_surface(f,:), Fs/(2*1000), rate_L_sample_freqs_kHz);
    %printf("f: %d %f %f %f\n", f, rate_K_sample_freqs_kHz(1), rate_L_sample_freqs_kHz(1), AmdB_(1));
    model_(f,1) = Wo; model_(f,2) = L; model_(f,3:(L+2)) = 10 .^ (AmdB_(f, 1:L)/20);
   end
endfunction


% PostFilter, has a big impact on speech quality after VQ.  When used
% on a mean removed rate K vector, it raises formants, and supresses
% anti-formants.  As it manipulates amplitudes, we normalise energy to
% prevent clipping or large level variations.  pf_gain of 1.2 to 1.5
% (dB) seems to work OK.  Good area for further investigations and
% improvements in speech quality.

function vec = post_filter(vec, sample_freq_kHz, pf_gain = 1.5, voicing)
    % vec is rate K vector describing spectrum of current frame
    % lets pre-emp before applying PF. 20dB/dec over 300Hz
    
    pre = 20*log10(sample_freq_kHz/0.3);
    vec += pre;

    levels_before_linear = 10 .^ (vec/20);
    e_before = sum(levels_before_linear .^2);

    vec *= pf_gain;
   
    levels_after_linear = 10 .^ (vec/20);
    e_after = sum(levels_after_linear .^2);
    gain = e_after/e_before;
    gaindB = 10*log10(gain);
    vec -= gaindB;

    vec -= pre;
endfunction


% construct energy quantiser table, and save to text file to include in C

function energy_q = create_energy_q
    energy_q = 10 + 40/16*(0:15);
endfunction

function save_energy_q(fn)
  energy_q = create_energy_q;
  f = fopen(fn, "wt");
  fprintf(f, "1 %d\n", length(energy_q));
  for n=1:length(energy_q)
    fprintf(f, "%f\n", energy_q(n));
  end
  fclose(f);
endfunction


% save's VQ in format that can be compiled by Codec 2 build system

function save_vq(vqset, filenameprefix)
  [Nvec order stages] = size(vqset);
  for s=1:stages
    fn = sprintf("%s_%d.txt", filenameprefix, s);
    f = fopen(fn, "wt");
    fprintf(f, "%d %d\n", order, Nvec);
    for n=1:Nvec
      for k=1:order
        fprintf(f, "% 8.4f ", vqset(n,k,s));
      end
      fprintf(f, "\n");
    end
    fclose(f);
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


function [diff_weighted weights error g min_ind] = search_vq_weighted(target, vq, weight_gain)
  [vq_rows vq_cols] = size(vq);

  weight_gain = 0.1; % I like this vairable name as it is funny

  % find mse for each vector

  error = g = zeros(1, vq_rows);
  diff = weights = diff_weighted = zeros(vq_rows, vq_cols);

  weights = max(0.1, weight_gain .* (target + 20));

  for i=1:vq_rows

    % work out gain for best match

    g(i) = sum((target - vq(i,:)).*weights)/vq_cols;

    % Find weighted difference.  This allocated more importance
    % (error) to samples with higher energy, and stops really low
    % level harmonics from having any impact.  Note addition in dB
    % is multiplication in linear

    diff(i,:) = target - vq(i,:) - g(i);
 
    diff_weighted(i,:) = diff(i,:) .* weights;

    % abs in dB is MSE in linear

    error(i) = mean(abs(diff_weighted(i,:)));
  end

  [mn min_ind] = min(error);
  
endfunction


% --------------------------------------------------------------------------------
% Experimental functions used for masking, piecewise models, not part of newamp1
% --------------------------------------------------------------------------------


function [maskdB_ maskdB_cyclic Dabs dk_ D1_ ind] = decimate_in_freq(maskdB, cyclic=1, k=7, vq)

    % Lets try to come up with a smoothed, cyclic model.  Replace
    % points from 3500 Hz to 4000Hz with a sequence that joins without
    % a step to points at the 0Hz end of the spectrum.  This will make
    % it more cyclical and make the DFT happier, less high freq
    % energy.  Yes, happier is an extremely technical term.

    L = length(maskdB);
    anchor = floor(7*L/8);
    xpts = [ anchor-1 anchor L+1 L+2];
    ypts = [ maskdB(anchor-1) maskdB(anchor) maskdB(1) maskdB(2)];
    mask_pp = splinefit(xpts, ypts, 1);
    maskdB_cyclic = [maskdB(1:anchor) ppval(mask_pp, anchor+1:L)];

    % Now DFT, truncating DFT coeffs to undersample

    if cyclic
      D = fft(maskdB_cyclic)/L;
    else
      D = fft(maskdB)/L;
    end
    Dabs = abs(D);                        % this returned for plotting

    % truncate D to rate k, convert to 2k length real vector for quantisation and transmission

    Dk = [0 D(2:k-1) real(D(k)) D(L-k+1:L)]; 
    dk = real(ifft(Dk));
    D1 = D(1);
       
    % quantisation

    if nargin == 4
       [res tmp vq_ind] = mbest(vq, dk, 4);
       D1_tab = 0:(60/15):60;
       assert(length(D1_tab) == 16);
       [tmp D1_ind] = quantise(D1_tab, D1);
       ind = [vq_ind D1_ind];
       [dk_ D1_] = index_to_params(ind, vq);
       %std(dk_ - dk)
    else
       dk_ = dk;
       D1_ = D1;
    end

    maskdB_ = params_to_mask(L, k, dk_, D1_);
endfunction



function [dk_ D1_] = index_to_params(ind, vq)
    [Nvec order stages] = size(vq);
    dk_ = zeros(1,order);
    for s=1:stages
      dk_ = dk_ + vq(ind(s),:,s);
    end
    D1_tab = 0:(60/15):60;
    D1_ = D1_tab(ind(stages+1));
endfunction


% decoder side

function maskdB_ = params_to_mask(L, k, dk_, D1_)

    anchor = floor(7*L/8);

    % convert quantised dk back to rate L magnitude spectrum

    Dk_ = fft(dk_);
    D_ = zeros(1,L);
    D_(1) = D1_;                      % energy seperately quantised
    D_(2:k-1) = Dk_(2:k-1);
    D_(L-k+1:L) = Dk_(k+1:2*k);
    d_ = L*ifft(D_);                  % back to spectrum at rate L
    maskdB_ = real(d_);
    
    % Finally fix up last 500Hz, taper down 10dB at 4000Hz

    xpts = [ anchor-1 anchor L];
    ypts = [ maskdB_(anchor-1) maskdB_(anchor) (maskdB_(anchor)-10)];
    mask_pp = splinefit(xpts, ypts, 1);
    maskdB_ = [maskdB_(1:anchor) ppval(mask_pp, anchor+1:L)];
endfunction



% determine cumulative mask, using amplitude of each harmonic.  Mask is
% sampled across L points in the linear domain

function maskdB = determine_mask(masker_amps_dB, masker_freqs_kHz, mask_sample_freqs_kHz, bark_model=1)

    % calculate and plot masking curve

    maskdB = -20*ones(1,length(mask_sample_freqs_kHz));
    for m=1:length(masker_freqs_kHz)
      maskdB = max(maskdB, schroeder(masker_freqs_kHz(m), mask_sample_freqs_kHz, bark_model) + masker_amps_dB(m)); 
      %maskdB = max(maskdB, parabolic_resonator(masker_freqs_kHz(m), mask_sample_freqs_kHz) + masker_amps_dB(m)); 
    end
end


% Sample mask as model for Am

function [maskdB Am_freqs_kHz] = mask_model(AmdB, Wo, L, bark_model=1)

    Am_freqs_kHz = (1:L)*Wo*4/pi;
    maskdB = determine_mask(AmdB, Am_freqs_kHz, Am_freqs_kHz, bark_model);
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

    %x1 = 0.5; x2 = 2;
    %y1 = 0.5; y2 = 1;
    x1 = 0.5; x2 = 3;
    y1 = 1; y2 = 3;

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


% Alternative mask function that has a gentler slope than schroeder.
% Idea is to get sharp formant definition, but also fill in gaps so we
% dont get chunks of spectrum coming and going

function [maskdB pp] = resonator(freq_tone_kHz, mask_sample_freqs_kHz)
  % note all frequencies in kHz

  f1  = 0.1; f2  = 3;
  bw1 = 0.1; bw2 = 0.1; 
  m = (bw2-bw1)/(log10(f2)-log10(f1));
  c = bw1 - m*log10(f1);

  Fs    = 8;
  slope = -12;     % filter falls off by this slope/octave

  maskdB = zeros(1, length(mask_sample_freqs_kHz));

  % frequency dependant bandwidth

  bw = m*log10(freq_tone_kHz) + c;
  printf("freq_tone_kHz: %f bw: %f\n", freq_tone_kHz, bw);

  % Design spline to set shape based on current bandwidth

  x = [-Fs/2 -bw/2 0 +bw/2 +Fs/2];
  delta = slope*log2(Fs/bw);                 % gain is delta down from -3dB to Fs/2
  y = [-3 + delta, -3, 0, -3, -3 + delta];
  pp = splinefit(x, y, 4);
  maskdB = ppval(pp, mask_sample_freqs_kHz - freq_tone_kHz);
endfunction


function maskdB = resonator_fast(freq_tone_kHz, mask_sample_freqs_kHz)

  % note all frequencies on kHz

  #{
  max_ind = length(pp_bw);
  ind = round(freq_tone_kHz/0.1);
  ind = min(ind, max_ind);
  ind = max(ind, 1);
  #}
  %printf("freq_tone_kHz: %f ind: %d\n", freq_tone_kHz, ind);
  [maskdB_res1 pp] = resonator(0.5, mask_sample_freqs_kHz);

  maskdB = ppval(pp, mask_sample_freqs_kHz - freq_tone_kHz);
  %maskdB = ppval(pp_bw(ind), mask_sample_freqs_kHz - freq_tone_kHz);
endfunction


% Alternative mask function that uses parabolas for fast computation

function maskdB = parabolic_resonator(freq_tone_kHz, mask_sample_freqs_kHz)

  % note all frequencies in kHz

  % bandwidth as a function of log(f)

  f1  = 0.5; f2  = 3;
  bw1 = 0.1; bw2 = 0.3; 
  m = (bw2-bw1)/(log10(f2)-log10(f1));
  c = bw1 - m*log10(f1);

  Fs = 8;
  slope = -18;

  % frequency dependant bandwidth

  if freq_tone_kHz < f1
    bw = bw1;
  else
    bw = m*log10(freq_tone_kHz) + c;
  end
  %printf("freq_tone_kHz: %f bw: %f\n", freq_tone_kHz, bw);

  % Design parabola to fit bandwidth

  a = -3/((bw/2)^2);
  %printf("freq_tone_kHz: %f bw: %f a: %f\n", freq_tone_kHz, bw, a);

  % Design straight line to fit slope
 
  delta = slope*log2(Fs/bw);                 % gain is delta down from -3dB to Fs/2
  m1 = 2*delta/Fs;

  maskdB_par  = a*((mask_sample_freqs_kHz - freq_tone_kHz).^2);
  maskdB_line = m1*abs(mask_sample_freqs_kHz - freq_tone_kHz) - 10;
  %indx = find(mask_sample_freqs_kHz < freq_tone_kHz);
  %maskdB_line(indx) = -50;

  maskdB = max(maskdB_par, maskdB_line);
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

function maskdB_ = decimate_frame_rate(model, decimate, f, frames)
    max_amp = 80;

    Wo = model(f,1);
    L = min([model(f,2) max_amp]);

    % determine frames that bracket the one we need to interp

    left_f = decimate*floor((f-1)/decimate)+1; 
    right_f = left_f + decimate;
    if right_f > frames
      right_f = left_f;
    end

    % determine fraction of each frame to use

    left_fraction  = 1 - mod((f-1),decimate)/decimate;
    right_fraction = 1 - left_fraction;

    printf("f: %d left_f: %d right_f: %d left_fraction: %3.2f right_fraction: %3.2f \n", f, left_f, right_f, left_fraction, right_fraction)

    % fit splines to left and right masks

    left_Wo = model(left_f,1);
    left_L = min([model(left_f,2) max_amp]);
    left_AmdB = 20*log10(model(left_f,3:(left_L+2)));
    left_sample_freqs_kHz = (1:left_L)*left_Wo*4/pi;

    right_Wo = model(right_f,1);
    right_L = min([model(right_f,2) max_amp]);
    right_AmdB = 20*log10(model(right_f,3:(right_L+2)));
    right_sample_freqs_kHz = (1:right_L)*right_Wo*4/pi;

    % determine mask for left and right frames, sampling at Wo for this frame

    sample_freqs_kHz = (1:L)*Wo*4/pi;
    maskdB_left = interp1(left_sample_freqs_kHz, left_AmdB, sample_freqs_kHz, "extrap");
    maskdB_right = interp1(right_sample_freqs_kHz, right_AmdB, sample_freqs_kHz, "extrap");

    maskdB_ = left_fraction*maskdB_left + right_fraction*maskdB_right;
endfunction


% plot some masking curves, used for working on masking filter changes

function plot_masking(bark_model=0);
  Fs = 8000;

  figure(1)
  mask_sample_freqs_kHz = 0.1:0.025:(Fs/1000)/2;
  %maskdB_s0 = schroeder(0.5, mask_sample_freqs_kHz, 0);
  %plot(mask_sample_freqs_kHz, maskdB_s0,';schroeder 0;');
  maskdB_s1 = schroeder(0.5, mask_sample_freqs_kHz, bark_model);
  plot(mask_sample_freqs_kHz, maskdB_s1,'g;schroeder 1;');
  #{
  maskdB_res = parabolic_resonator(0.5, mask_sample_freqs_kHz);
  plot(mask_sample_freqs_kHz, maskdB_res,'r;resonator;');
  #}
  hold on;

  for f=0.5:0.5:3
    %maskdB_s0 = schroeder(f, mask_sample_freqs_kHz, 0);
    %plot(mask_sample_freqs_kHz, maskdB_s0);
    maskdB_s1 = schroeder(f, mask_sample_freqs_kHz, bark_model);
    plot(mask_sample_freqs_kHz, maskdB_s1,'g');
    #{
    maskdB_res = parabolic_resonator(f, mask_sample_freqs_kHz);
    plot(mask_sample_freqs_kHz, maskdB_res,'r;resonator;');
    #}
  end
  hold off;

  axis([0.1 4 -30 0])
  grid

  #{
  %pp_bw = gen_pp_bw;
  figure(2)
  clf;
  maskdB_res = resonator(0.5, mask_sample_freqs_kHz);
  plot(mask_sample_freqs_kHz, maskdB_res);  
  hold on;
  maskdB_res_fast = resonator_fast(0.5, mask_sample_freqs_kHz);
  plot(mask_sample_freqs_kHz, maskdB_res_fast, "g");  
  maskdB_par = parabolic_resonator(0.5, mask_sample_freqs_kHz);
  plot(mask_sample_freqs_kHz, maskdB_par, "r");  
  hold off;
  axis([0 4 -80 0]);
  grid
  #}
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

function D1_log = decode_from_bit_stream(samname, ber = 0, bit_error_mask = ones(28,1))
  max_amp = 80;
  bits_per_param = [6 1 8 8 4 1];
  assert(sum(bits_per_param) == 28);
  load vq;
  k = 10;
  dec_in_time = 1;
  train = 0;
  decimate = 4;
  synth_phase = 1;

  Am_out_name = sprintf("%s_am.out", samname);
  Aw_out_name = sprintf("%s_aw.out", samname);
  bit_stream_name = strcat(samname,".bit");
  faw = fopen(Aw_out_name,"wb"); 
  fam  = fopen(Am_out_name,"wb"); 
  faw = fopen(Aw_out_name,"wb"); 

  Wo_out_name = sprintf("%s_Wo.out", samname);
  fWo  = fopen(Wo_out_name,"wb"); 

  % read in bit stream and convert to ind_log[]

  ind_log = [];
  fbit  = fopen(bit_stream_name, "rb"); 
  bits_per_frame = sum(bits_per_param);
  nind = length(bits_per_param);
  nerr = 0; nbits = 0;
  [frame nread] = fread(fbit, sum(bits_per_param), "uchar");
  while (nread == bits_per_frame)

    % optionally introduce bit errors

    error_bits = rand(sum(bits_per_param), 1) < ber;
    error_bits_masked = bitand(error_bits, bit_error_mask);
    frame = bitxor(frame, error_bits_masked);
    nerr += sum(error_bits_masked);
    nbits += sum(bits_per_param);

    % read a frame, convert to indexes

    nbit = 1;
    ind = [];
    for i=1:nind
      field = frame(nbit:nbit+bits_per_param(i)-1);
      nbit += bits_per_param(i);
      ind = [ind bits_to_index(field, bits_per_param(i))];
    end
    ind_log = [ind_log; ind];
    [frame nread] = fread(fbit, sum(bits_per_param), "uchar");
  endwhile
  fclose(fbit);
  printf("nerr: %d nbits: %d %f\n", nerr, nbits, nerr/nbits);
  
  % convert ind_log to modem params

  frames = 4*length(ind_log);
  model_ = zeros(frames, max_amp+2);
  v      = zeros(frames,1);
  D1_log = [];

  fdec = 1;
  for f=1:4:frames
    ind_Wo = ind_log(fdec,1);

    Wo = decode_log_Wo(ind_Wo, 6);
    L = floor(pi/Wo);
    L = min([L max_amp-1]);
    model_(f,1) = Wo;
    model_(f,2) = L;
    
    v1 = ind_log(fdec,2); 
    if (fdec+1) < length(ind_log)
      v5 = ind_log(fdec+1,2);
    else
      v5 = 0;
    end
    v(f:f+3) = est_voicing_bits(v1, v5);

    ind_vq = ind_log(fdec,3:5) + 1;
    [dk_ D1_] = index_to_params(ind_vq, vq);
    D1_log = [D1_log; D1_];
    maskdB_ = params_to_mask(L, k, dk_, D1_);
    Am_ = zeros(1,max_amp);
    Am_ = 10 .^ (maskdB_(1:L)/20); 
    model_(f,3:(L+2)) = Am_;

    fdec += 1;
  end

  % decoder loop -----------------------------------------------------

  if train
    % short circuit decoder
    frames = 0;
  end

  % run post filter ahead of time so dec in time has post filtered frames to work with

  for f=1:frames
    model_(f,:) = post_filter(model_(f,:));
  end

  for f=1:frames
    %printf("frame: %d\n", f);
    L = min([model_(f,2) max_amp-1]);
    Wo = model_(f,1);
    Am_ = model_(f,3:(L+2));

    maskdB_ = 20*log10(Am_);

    if dec_in_time
      % decimate mask samples in time

      [maskdB_ Wo L] = decimate_frame_rate2(model_, decimate, f, frames);
      model_(f,1) = Wo;
      model_(f,2) = L;
    end

    Am_ = zeros(1,max_amp);
    Am_(2:L) = 10 .^ (maskdB_(1:L-1)/20);  % C array doesnt use A[0]
    fwrite(fam, Am_, "float32");
    fwrite(fWo, Wo, "float32");

    if synth_phase

      % synthesis phase spectra from magnitiude spectra using minimum phase techniques

      fft_enc = 512;
      model_(f,3:(L+2)) = 10 .^ (maskdB_(1:L)/20);
      phase = determine_phase(model_, f);
      assert(length(phase) == fft_enc);
      Aw = zeros(1, fft_enc*2); 
      Aw(1:2:fft_enc*2) = cos(phase);
      Aw(2:2:fft_enc*2) = -sin(phase);
      fwrite(faw, Aw, "float32");    
    end
  end

  fclose(fam);
  fclose(fWo);
  if synth_phase
    fclose(faw);
  end

  % save voicing file
  
  v_out_name = sprintf("%s_v.txt", samname);
  fv  = fopen(v_out_name,"wt"); 
  for f=1:length(v)
    fprintf(fv,"%d\n", v(f));
  end
  fclose(fv);

endfunction



% decimate frame rate of mask, use linear interpolation in the log domain 

function [maskdB_ Wo L] = decimate_frame_rate2(model, decimate, f, frames)
    max_amp = 80;

    % determine frames that bracket the one we need to interp

    left_f = decimate*floor((f-1)/decimate)+1; 
    right_f = left_f + decimate;
    if right_f > frames
      right_f = left_f;
    end

    % determine fraction of each frame to use

    left_fraction  = 1 - mod((f-1),decimate)/decimate;
    right_fraction = 1 - left_fraction;

    % printf("f: %d left_f: %d right_f: %d left_fraction: %3.2f right_fraction: %3.2f \n", f, left_f, right_f, left_fraction, right_fraction)

    % fit splines to left and right masks

    left_Wo = model(left_f,1);
    left_L = min([model(left_f,2) max_amp]);
    left_AmdB = 20*log10(model(left_f,3:(left_L+2)));
    left_mask_sample_freqs_kHz = (1:left_L)*left_Wo*4/pi;
    
    right_Wo = model(right_f,1);
    right_L = min([model(right_f,2) max_amp]);
    right_AmdB = 20*log10(model(right_f,3:(right_L+2)));
    right_mask_sample_freqs_kHz = (1:right_L)*right_Wo*4/pi;

    % printf("  right_Wo: %f left_Wo: %f  right_L: %d  left_L %d\n",right_Wo,left_Wo,right_L,left_L);
    printf("%f %f\n", left_AmdB(left_L), right_AmdB(right_L));

    maskdB_left_pp = splinefit(left_mask_sample_freqs_kHz, left_AmdB, left_L);
    maskdB_right_pp = splinefit(right_mask_sample_freqs_kHz, right_AmdB, right_L);

    % determine mask for left and right frames, sampling at Wo for this frame

    Wo = left_fraction*left_Wo + right_fraction*right_Wo;
    L = floor(pi/Wo);
    %Wo = model(f,1); L = model(f,2);

    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    maskdB_left = ppval(maskdB_left_pp, mask_sample_freqs_kHz);
    maskdB_right = ppval(maskdB_right_pp, mask_sample_freqs_kHz);

    maskdB_ = left_fraction*maskdB_left + right_fraction*maskdB_right;
endfunction

#{
function amodel = post_filter(amodel)
    max_amp = 80;

    % post filter 

    L = min([amodel(2) max_amp-1]);
    Wo = amodel(1);
    Am_ = amodel(3:(L+2));
    AmdB_ = 20*log10(Am_);
    AmdB_pf = AmdB_*1.5;
    AmdB_pf += max(AmdB_) - max(AmdB_pf);
    amodel(3:(L+2)) = 10 .^ (AmdB_pf(1:L)/20);
endfunction
#}


% Given a matrix with indexes on each row, convert to a bit stream and
% write to file.  We only write every 4th frame due to DIT

function write_bit_stream_file(fn, ind_log, bits_per_param)
  fbit  = fopen(fn,"wb"); 
  decimate = 4;

  % take a row of quantiser indexes, convert to bits, save to file

  [frames nind] = size(ind_log);
  for f=1:decimate:frames
    frame_of_bits = [];
    arow = ind_log(f,:);
    for i=1:nind
      %printf("i: %d bits_per_param: %d\n", i, bits_per_param(i));
      some_bits = index_to_bits(arow(i), bits_per_param(i));
      frame_of_bits = [frame_of_bits some_bits];
    end
    fwrite(fbit, frame_of_bits, "uchar");
  end
  fclose(fbit);
endfunction



% determine 4 voicing bits based on 2 decimated voicing bits

function [v] = est_voicing_bits(v1, v5)
  if v1 == v5
    v(1:4) = v1;
  else
    v(1:2) = v1;
    v(3:4) = v5;
  end
endfunction


function [AmdB_ residual fvec fvec_ amps] = piecewise_model(AmdB, Wo, vq, vq_m) 
    L = length(AmdB);
    l1000 = floor(L/4);     
    AmdB_ = ones(1,L);
    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;

    % fit a resonator to max of first 300 - 1000 Hz

    fmin = 0.150;
    lmin = floor(L*fmin/4);
    [mx mx_ind] = max(AmdB(lmin+1:l1000));
    amp(1) = mx;
    mx_ind += lmin;
    AmdB_ = parabolic_resonator(mx_ind*Wo*4/pi, mask_sample_freqs_kHz) + mx;  
    fr1 = mx_ind*Wo*4/pi;

    % fit a 2nd resonator, must be above 1000Hz

    fmin = 1;
    lmin = round(fmin*L/4);

    [mx mx_ind] = max(AmdB(lmin+1:L));
    amp(2) = mx;
    mx_ind += lmin;
    AmdB_ = max(AmdB_, parabolic_resonator(mx_ind*Wo*4/pi, mask_sample_freqs_kHz) + mx);  
    fr2 = mx_ind*Wo*4/pi;

    % fit a third resonator, must be +/- 300 Hz after 2nd resonator

    residual = AmdB - AmdB_;
    keep_out = [1:lmin];
    lmax = round(L*3500/4000);
    keep_out = [1:lmin lmax:L];
    residual(keep_out) = -40;

    fr2 = mx_ind*Wo*4/pi;
    fmin = fr2 - 0.300;
    fmax = fr2 + 0.300;   
    lmin = max(1, round(L*fmin/4));
    lmax = min(L, round(L*fmax/4));
    keep_out = [keep_out lmin:lmax];

    residual = AmdB;
    residual(keep_out) = -40;

    if 0
    figure(3); clf;
    subplot(211)
    plot(mask_sample_freqs_kHz, residual);
    end

    [mx mx_ind] = max(residual);
    amp(3) = AmdB(mx_ind);
    AmdB_ = max(AmdB_, parabolic_resonator(mx_ind*Wo*4/pi, mask_sample_freqs_kHz) + amp(3));  
    fr3 = mx_ind*Wo*4/pi;

    % 4th resonator 

    fmin = fr3 - 0.300;
    fmax = fr3 + 0.300;
    
    lmin = max(1, round(L*fmin/4));
    lmax = min(L, round(L*fmax/4));
    keep_out = [keep_out lmin:lmax];

    residual = AmdB - AmdB_;
    residual(keep_out) = -40;

    [mx mx_ind] = max(residual);
    amp(4) = AmdB(mx_ind);
    AmdB_ = max(AmdB_, parabolic_resonator(mx_ind*Wo*4/pi, mask_sample_freqs_kHz) + amp(4));  
    fr4 = mx_ind*Wo*4/pi;

    if 0
    subplot(212)
    plot(mask_sample_freqs_kHz, residual);
    end

    printf("\nfr1: %f fr2: %f fr3: %f fr4: %f\n", fr1, fr2, fr3, fr4);
    [fvec fvec_ind] = sort([fr1 fr2 fr3 fr4]);
    amps = amp(fvec_ind(1:4));

    fvec_ = zeros(1, 4);

    #{
    % optional VQ of frequencies
 
    if nargin == 4
      AmdB_ = ones(1,L);
      [mes fvec_ ind] = mbest(vq, fvec, vq_m);
      for i=1:4
        an_amp = amp(fvec_ind(i));
        AmdB_ = max(AmdB_, parabolic_resonator(fvec_(i), mask_sample_freqs_kHz) + an_amp);
      end
    end
    #}

    % optional VQ of amplitudes
 
    if nargin == 4
      AmdB_ = ones(1,L);
      %amps_(1) = amps(1);
      %[mes tmp ind] = mbest(vq,  amps(2:4) -  amps_(1), vq_m);
      %amps_(2:4) = amps_(1) + tmp;
      [mes amps_ ind] = mbest(vq,  amps, vq_m);
      amps-amps_
      for i=1:4
        AmdB_ = max(AmdB_, parabolic_resonator(fvec(i), mask_sample_freqs_kHz) + amps_(i));
      end
    end

    %amps = amps(2:4) - amps(1);
endfunction


% find best place for resonator by closed loop min MSE search

function lmin = abys(AmdB_, AmdB, Wo, L, mask_sample_freqs_kHz)
  lstart = round(L/4);
  lmin = lstart;
  emin = 1E6;

  printf("lstart: %d L: %d\n", lstart, L);

  figure(3);
  subplot(211)
  plot(mask_sample_freqs_kHz*1000, AmdB,'r+-');
  
  e = zeros(1,L);
  for l=lstart:L

    % calc mse

    f_l = l*Wo*4/pi;
    AmdB_l = max(AmdB_, parabolic_resonator(f_l, mask_sample_freqs_kHz) + AmdB(l));
    hold on;
    if l == 23
      plot(mask_sample_freqs_kHz*1000, AmdB_l,'c');
    end
    hold off;
    e(l) = sum((AmdB_l - AmdB) .^ 2);
    %printf("l: %5d f_l: %4.3f e: %4.0f emin: %4.0f lmin: %5d\n", l, f_l, emin, lmin);
    printf("l: %5d f_l: %4.3f e: %4.0f emin: %4.0f lmin: %5d\n", l, f_l, e(l), emin, lmin);
    if e(l) < emin
       emin = e(l);
       lmin = l;
    end
  end

  subplot(212)
  plot(mask_sample_freqs_kHz*1000, e)
endfunction


function rate_K_surface_no_slope = remove_slope(rate_K_surface)
  [frames K] = size(rate_K_surface);
  rate_K_surface_no_slope = zeros(frames,K);
  for f=1:frames
    [gradient intercept] = linreg(1:K, rate_K_surface(f,:), K);
    printf("f: %d gradient: %f intercept: %f\n", f, gradient, intercept);
    rate_K_surface_no_slope(f,:) = rate_K_surface(f,:) - (intercept + gradient*(1:K));
  end
endfunction
