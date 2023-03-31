% newamp_700c.m
%
% Copyright David Rowe 2017
% This program is distributed under the terms of the GNU General Public License
% Version 2
%
% Library of Octave functions for rate K, mel spaced
% vector quantisation of spectral magnitudes used in Codec 2 700C mode.

1;
melvq; % mbest VQ functions

% --------------------------------------------------------------------------------
% Functions used by rate K mel work
% --------------------------------------------------------------------------------

% General 2nd order parabolic interpolator.  Used splines orginally,
% but this is much simpler and we don't need much accuracy.  Given two
% vectors of points xp and yp, find interpolated values y at points x

function y = interp_para(xp, yp, x)
  assert( (length(xp) >=3) && (length(yp) >= 3) );

  y = zeros(1,length(x));
  k = 1;
  for i=1:length(x)
    xi = x(i);

    % k is index into xp of where we start 3 points used to form parabola

    while ((xp(k+1) < xi) && (k < (length(xp)-2)))
      k++;
    end

    x1 = xp(k); y1 = yp(k); x2 = xp(k+1); y2 = yp(k+1); x3 = xp(k+2); y3 = yp(k+2);
    %printf("k: %d i: %d xi: %f x1: %f y1: %f\n", k, i, xi, x1, y1);

    a = ((y3-y2)/(x3-x2)-(y2-y1)/(x2-x1))/(x3-x1);
    b = ((y3-y2)/(x3-x2)*(x2-x1)+(y2-y1)/(x2-x1)*(x3-x2))/(x3-x1);

    y(i) = a*(xi-x2)^2 + b*(xi-x2) + y2;
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

  Gdbfk = interp_para(rate_L_sample_freqs_kHz, AmdB, sample_freqs_kHz);

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


function rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K)
  mel_start = ftomel(200); mel_end = ftomel(3700);
  step = (mel_end-mel_start)/(K-1);
  mel = mel_start:step:mel_end;
  rate_K_sample_freqs_Hz = 700*((10 .^ (mel/2595)) - 1);
  rate_K_sample_freqs_kHz = rate_K_sample_freqs_Hz/1000;
endfunction


function [rate_K_surface rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model, K)
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K);
  rate_K_surface = resample_const_rate_f(model, rate_K_sample_freqs_kHz);
endfunction


function rate_K_sample_freqs_kHz = lin_sample_freqs_kHz(K)
  lin_start = 100; lin_end = 3700;
  step = (lin_end-lin_start)/(K-1);
  rate_K_sample_freqs_Hz = lin_start:step:lin_end;
  rate_K_sample_freqs_kHz = rate_K_sample_freqs_Hz/1000;
endfunction


function [rate_K_surface rate_K_sample_freqs_kHz] = resample_const_rate_f_lin(model, K, resampler='spline')
  rate_K_sample_freqs_kHz = lin_sample_freqs_kHz(K);
  rate_K_surface = resample_const_rate_f(model, rate_K_sample_freqs_kHz, clip_en = 0, resampler);
endfunction


% Resample Am from time-varying rate L=floor(pi/Wo) to fixed rate K.  This can be viewed
% as a 3D surface with time, freq, and ampitude axis.

function rate_K_surface = resample_const_rate_f(model,
                                                rate_K_sample_freqs_kHz,
                                                clip_en=1,
                                                resampler='para')
  % convert rate L=pi/Wo amplitude samples to fixed rate K

  max_amp = 80;
  [frames col] = size(model);
  K = length(rate_K_sample_freqs_kHz);
  rate_K_surface = zeros(frames, K);

  for f=1:frames
    Wo = model(f,1);
    L = min([model(f,2) max_amp-1]);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);

    % clip between peak and peak -50dB, to reduce dynamic range

    if clip_en
      AmdB_peak = max(AmdB);
      AmdB(find(AmdB < (AmdB_peak-50))) = AmdB_peak-50;
    end

    rate_L_sample_freqs_kHz = (1:L)*Wo*4/pi;

    if strcmp(resampler, 'para')
      rate_K_surface(f,:) = interp_para(rate_L_sample_freqs_kHz, AmdB, rate_K_sample_freqs_kHz);
    else
      rate_K_surface(f,:) = interp1(rate_L_sample_freqs_kHz, AmdB, rate_K_sample_freqs_kHz, "spline", "extrap");
    end
    %printf("\r%d/%d", f, frames);
  end
  %printf("\n");
endfunction


% Take a rate K surface and convert back to time varying rate L

function [model_ AmdB_] = resample_rate_L(model, rate_K_surface, rate_K_sample_freqs_kHz, resampler='para')
  max_amp = 80;
  [frames col] = size(model);

  model_ = zeros(frames, max_amp+2);
  for f=1:frames
    Wo = model(f,1);
    L = model(f,2);
    rate_L_sample_freqs_kHz = (1:L)*Wo*4/pi;

    % back down to rate L

    if strcmp(resampler, 'para')
      AmdB_ = interp_para([ 0 rate_K_sample_freqs_kHz 4], [0 rate_K_surface(f,:) 0], rate_L_sample_freqs_kHz);
    else
      AmdB_ = interp1([0 rate_K_sample_freqs_kHz 4], [0 rate_K_surface(f,:) 0], rate_L_sample_freqs_kHz, "spline", 0);
    end
    model_(f,1) = Wo; model_(f,2) = L; model_(f,3:(L+2)) = 10 .^ (AmdB_(1:L)/20);
   end
endfunction


% Post Filter, has a big impact on speech quality after VQ.  When used
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


% Equaliser in front of EQ, see vq_700c_eq.m for development version

function [rate_K_vec eq] = front_eq(rate_K_vec, eq)
  [tmp K] = size(rate_K_vec);
  ideal = [ 8 10 12 14 14*ones(1,K-1-4) -20];
  gain = 0.02;
  update = rate_K_vec - ideal;
  eq = (1-gain)*eq + gain*update;
  eq(find(eq < 0)) = 0;
endfunction

% ----------------------------------------------------------------------
% Functions added in "Rate K" experimental R&D campaign from August 2022
% ----------------------------------------------------------------------

function fHz = warp(k, K)
  mel_start = ftomel(200); mel_end = ftomel(3700);
  step = (mel_end-mel_start)/(K-1);
  mel = mel_start + (k-1)*step;
  fHz = 700*((10 .^ (mel/2595)) - 1);
endfunction


function k = warp_inv(fHz, K)
  mel_start = ftomel(200); mel_end = ftomel(3700);
  step = (mel_end-mel_start)/(K-1);
  mel = ftomel(fHz) - step;
  k = (ftomel(fHz) - mel_start)/step + 1;
endfunction


function h = generate_filter(m,f0,L,Nb)
  g = zeros(1,L);
  fb = m*f0;
  b = warp_inv(fb,Nb);
  st = round(warp(b-1,Nb)/f0); st = max(st,1);
  en = round(warp(b+1,Nb)/f0); en = min(en,L);
  %printf("fb: %f b: %f warp(b-1): %f warp(b+1): %f st: %d en: %d\n", fb, b, warp(b-1,Nb), warp(b+1,Nb), st, en);
  for k=st:m-1
    g(k) = (k-st)/(m-st);
  end
  g(m) = 1;
  for k=m+1:en
    g(k) = (en-k)/(en-m);
  end
  g_sum = sum(g);
  h = g/g_sum;
endfunction

% make sure nothing breaks when we generate a bunch of filters
function test_filters
  Fs=8000; f0=200; L=floor(Fs/(2*f0));
  for m=1:L
    generate_filter(m,f0,L,20);
  end
end

% Given a magnitude spectrum, compute the LPC coeffcients
% w and mag row vectors
function [a G] = lpc_from_mag_spectrum(w, mag, order)
  assert(isrow(w));
  assert(isrow(mag));
  assert(length(w)==length(mag));
  M = length(w);
  P = mag.^2;
  R = zeros(order+1,1); % R(i) actually stored in R(i+1) because Octave
  for i=0:order
    R(i+1) = (1/M)*P*cos(i*w)';
  end
  R_toeplitz = toeplitz(R(1:order));
  a = -inv(R_toeplitz)*R(2:order+1);
  G = sqrt(R(1) + a'*R(2:order+1));
  a = [1 a'];
end

function test_lpc_from_mag_spectrum
  step=pi/100;
  w = 0:step:pi-step;
  omega = pi/4; beta=0.9;
  a_target =  [1 -2*beta*cos(omega) beta*beta];
  mag = abs(freqz(1, a_target, w));
  [a_est G] = lpc_from_mag_spectrum(w, mag, 2);
  mag_est = abs(freqz(1, a_est, w));
  SD = mean((20*log10(mag)-20*log10(mag_est)).^2);
  assert(SD < 0.1);
  assert(abs(G-1) < 0.1);
  clf; plot(w,mag); hold on; plot(w,mag_est); hold off;
end

% Synthesised phase0 model using Hilbert Transform
function phase0 = synth_phase_from_mag(mag_sample_freqs_kHz, mag_dB, Fs, Wo, L, postfilter_en)
  Nfft=512;
  sample_freqs_kHz = (Fs/1000)*[0:Nfft/2]/Nfft;  % fft frequency grid (nonneg freqs)
  Gdbfk = interp1([0 mag_sample_freqs_kHz 4], [0 mag_dB 0], sample_freqs_kHz, "spline", "extrap");
  if postfilter_en
    % optional "all pass" or "phase" postfiltering
    % in a minimum phase system, group delay is proportional to slope of mag spectrum
    Gdbfk *= 1.5;
  end
  [phase_ht s] = mag_to_phase(Gdbfk, Nfft);
  phase0 = zeros(1,L);
  for m=1:L
    b = round(m*Wo*Nfft/(2*pi));
    phase0(m) = phase_ht(b);
  end
end

function [YdB SdB] = amplitude_postfilter(rate_Lhigh_sample_freqs_kHz, YdB, Fs, F0high, eq=0)
  % straight line fit to YdB to estimate spectral slope SdB
  w = 2*pi*rate_Lhigh_sample_freqs_kHz*1000/Fs;
  st = max(round(200/F0high),1);
  en = min(round(3700/F0high), length(w));
  %printf("st: %d en: %d\n", st, en);
  [m b] = linreg(w(st:en),YdB(st:en),en-st+1);
  SdB = w*m+b;

  Y = 10 .^ (YdB/20);
  Y_energy1 = sum(Y .^ 2);

  % remove slope and expand dynamic range
  YdB -= SdB;
  m1 = (2.0-1.2)/pi; c1 = 1.2;
  YdB = YdB .* (m1*w+c1);
  if eq == 0, YdB += SdB; end
  if eq == 1, YdB += 0.5*SdB; end
  if eq == 2
    %dynamic_range = 30;
    %mx = max(YdB); YdB = max(YdB, mx-dynamic_range);
    %printf("  b: %f m: %f\n",b, m);
    
    % whiten voiced speech (voiced speech has -6dB/octave fall off)
    if (m > 0)
      YdB += w*m;
    end
  end
  
  % normalise energy
  Y = 10 .^ (YdB/20);
  Y_energy2 = sum(Y .^ 2);
  YdB += 10*log10(Y_energy1/Y_energy2);
end

function str = papr(s)
  papr_dB = 10*log10(max(abs(s).^2)/mean(abs(s).^2));
  str = sprintf("CPAPR: %3.1f dB", papr_dB);
end

function [Eq best_i gmin] = weighted_search(vq_stage1, target);
  mx = max(target);
  w = (0.75/30)*(target-mx) + 1.0;
  w2 = w .^ 2;
  I = length(vq_stage1);
  Emin = 1E32;
  for i=1:I
    g = sum((target-vq_stage1(i,:)).*w2)/sum(w2);
    E = sum( ((target - vq_stage1(i,:) - g) .* w) .^ 2 );
    if E < Emin
      best_i = i;
      Emin = E;
      gmin = g;
    end
  end
  Eq = Emin/length(target);
end

function f = eq_cand_b(f, target, vq, delta=0.01)
  g = mean(target);
  [res target_ ind] = mbest(vq, target-f-g, 1);
  f += 2*delta*(target-f-g-target_);
end

% normalises the energy in AmdB_rate2 to be the same as AmdB_rate1
function AmdB_rate2_hat = norm_energy(AmdB_rate1, AmdB_rate2)
  c = sum(10 .^ (AmdB_rate1/10))/sum(10 .^ (AmdB_rate2/10));
  AmdB_rate2_hat = AmdB_rate2 + 10*log10(c);
end

% combine interpolation with energy normalisation
function YIdB_norm = interp1_norm(XkHz,YdB,XIkHz)
  YIdB = interp1([0 XkHz 4],[0 YdB 0],XIkHz,"spline", "extrap");
  YIdB_norm = norm_energy(YdB,YIdB);
end

function test_norm_energy
  Fs = 8000; L = 80; F0 = (Fs/2)/L; K=20;
  Am_freqs_kHz = (1:L-1)/1000;
  AmdB = ones(1,L-1);
 
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K);
  AmdB_rate_K_hat = interp1_norm(Am_freqs_kHz, AmdB,
                    rate_K_sample_freqs_kHz);
  E_L = sum(10 .^ (AmdB/10));
  E_K = sum(10 .^ (AmdB_rate_K_hat/10));
  if abs(E_L-E_K) < 1E-3
    printf("PASS\n");
  else
    printf("FAIL\n");
  end
end

% Returns mu-law compressed frame energy. F_dB is scaled to 0dB when E_db = E_max_dB
function F_dB = mulaw_comp(E_dB)
  E_lin = 10 .^ (E_dB/10);
  E_max_dB = 60;
  s_max = 10 .^ (E_max_dB/20);
  x = sqrt(E_lin)/s_max;
  mu = 255;
  F = log(1+mu*x)/log(1+mu);
  F_dB = 20*log10(F);
endfunction

function E_dB = mulaw_decomp(F_dB)
  y = 10 .^ (F_dB/20);
  mu = 255;
  F_hat = (((1+mu).^y) - 1)/mu;
  E_max_dB = 60;
  s_max = 10 .^ (E_max_dB/20);
  E_hat = F_hat*s_max;
  E_dB = 20*log10(E_hat);
endfunction

% F and F_hat should be exact inverses when not quantised 
function test_mulaw
  E_dB=-20:80; 
  F_dB = mulaw_comp(E_dB);
  E_hat_dB = mulaw_decomp(F_dB); 
  figure(1); clf; plot(E_dB,F_dB); grid;
  xlabel('E (dB)'); ylabel('F (dB)');
  figure(2); clf; plot(E_dB,E_hat_dB); grid;
  xlabel('E (dB)'); ylabel('E\_hat (dB)');
  hold on; plot(E_dB, E_dB - E_hat_dB,'r-'); hold off;
  assert(abs(E_dB - E_hat_dB) < 0.001);
endfunction

function test_mulaw_quant
  E_dB=-20:0.1:80; 
  F_dB = mulaw_comp(E_dB);
  levels = -20 + (20/16)*(0:15)
  Q_dB = zeros(1,length(E_dB));
  for i=1:length(E_dB)
    Q_dB(i) = quantise(levels, F_dB(i));
  end
  E_hat_dB = mulaw_decomp(Q_dB); 
  figure(1); clf; plot(E_dB,F_dB); grid;
  hold on; plot(E_dB,Q_dB,'r-'); hold off;
  xlabel('E (dB)'); ylabel('F (dB)');
  figure(2); clf; plot(E_dB,E_hat_dB); grid;
  xlabel('E (dB)'); ylabel('E\_hat (dB)');
  assert(abs(E_dB - E_hat_dB) < 0.001);
endfunction

function vec_out = mean_compression(vec_in, noise_gate, thresh, slope, maximum, makeup_gain)
  x = mean(vec_in,2);
  y = zeros(length(x),1);
  for i=1:length(x)
    if (x(i) > noise_gate) && (x(i) < thresh)
      y(i) = x(i) + makeup_gain;;
    elseif (x(i) >= thresh) && (x(i) < maximum)
      y(i) = thresh + (x(i)-thresh)*slope + makeup_gain;
    elseif x(i) >= maximum
      y(i) = thresh + (maximum-thresh)*slope + makeup_gain;
    end
   end
  vec_out = vec_in - x + y;
end

function test_mean_compression
  mean_in = -20:80;
  vec_in = ones(length(mean_in),K=20);
  vec_in = diag(mean_in)*vec_in;
  vec_out = mean_compression(vec_in, 0, 20, 0.5, 60, 60-40);
  figure(1); clf; plot(vec_in, vec_out); grid;
end

% energy compression using two stage piecewise linear curve
function y = piecewise_compressor(x1,y1,x2,y2,x3,y3,x)
  m1 = (y2-y1)/(x2-x1); c1 = y1 - m1*x1;
  m2 = (y3-y2)/(x3-x2); c2 = y2 - m2*x2;
  y = zeros(length(x),1);
  for i =1:length(x)
    if x(i) < x1
      y(i) = 0;
    elseif x(i) < x2
      y(i) = m1*x(i) + c1;
    elseif x(i) < x3
      y(i) = m2*x(i) + c2;
    else
      y(i) = y3;
    end  
  end   
endfunction

function test_piecewise_compressor
  x=-20:80; 
  y = piecewise_compressor(20,36,60,76,80,76,x);
  figure(1); clf; plot(x,y); grid;
endfunction

% some plot to explore relationship between energy and mean
function test_energy_and_mean(B)
  mean_dB = mean(B,2);
  E_dB = 10*log10(sum(10 .^ (B(:,:)/10),2));
  figure(1); clf;
  subplot(211); hist(mean_dB,20); title('mean (dB)');
  subplot(212); hist(E_dB,20); title('Energy (dB)');
  figure(2); clf; plot(mean_dB, E_dB,'+'); grid;
  xlabel('mean (dB)'); ylabel('E (dB)');
  figure(3);
  mx = max(mean_dB); Nsteps=25;
  cdf = empirical_cdf(mx*(1:Nsteps)/Nsteps,mean_dB);
  plot(mx*(1:Nsteps)/Nsteps,cdf); title('CDF Estimates'); grid;
endfunction

