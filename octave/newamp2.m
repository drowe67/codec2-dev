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
#}

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

    model_(f,1) = Wo; model_(f,2) = L; model_(f,3:(L+2)) = 10 .^ (AmdB_(f, 1:L)/20);
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


% Given a vector of rate K samples, huffman encodes/decodes delta
% amplitude, returning quantised samples

function rate_K_vec_ = huffman_quantise_rate_K(rate_K_vec)
  K = length(rate_K_vec);

  % whole thing is quantised to 6dB steps, as that doesn't seem to
  % introduce much distortion
  
  rate_K_vec = 6*round(rate_K_vec/6);

  % start with k=3, around 250Hz, we assume that's quantised as the
  % mean frame energy, as samples before that might be stuck in the
  % HPF.

  rate_K_vec_no_mean = rate_K_vec - rate_K_vec(3);
  rate_K_vec_no_mean_ = zeros(1,K);

  % huffman encoding od differences, ignoring m=3, using the following table
  %
  % 00    0
  % 10   +6
  % 11   -6
  % 010 -12
  % 011 +12

  levels = [0 6 -6 -12 12]; symbols = {[0 0],[1 0],[1 1],[0 1 0;],[0 1 1]};

  % move backwards to get target for first two samples
  
  bits = [];
  [quant_out best_i] = quantise(levels, rate_K_vec_no_mean(2));
  bits = [bits symbols{best_i}];
  rate_K_vec_no_mean_(2) = quant_out;
  
  [quant_out best_i] = quantise(levels, rate_K_vec_no_mean(1) - rate_K_vec_no_mean(2));
  bits = [bits symbols{best_i}];
  rate_K_vec_no_mean_(1) = quant_out + rate_K_vec_no_mean_(2);

  % then forwards for rest of the target samples
  
  for m=4:K
    target = rate_K_vec_no_mean(m) - rate_K_vec_no_mean_(m-1);
    [quant_out best_i] = quantise(levels, target);
    bits = [bits symbols{best_i}];
    rate_K_vec_no_mean_(m) = quant_out + rate_K_vec_no_mean_(m-1);    
  end

  printf("%d bits\n", length(bits));
  rate_K_vec_ = rate_K_vec_no_mean_ + rate_K_vec(3);
endfunction
