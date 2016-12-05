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

1;
melvq; % mbest VQ functions


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


% quantise input sample to nearest value in table, optionally return bianry code

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


% determine cumulative mask, using amplitude of each harmonic.  Mask is
% sampled across L points in the linear domain

function maskdB = determine_mask(masker_amps_dB, masker_freqs_kHz, mask_sample_freqs_kHz)

    % calculate and plot masking curve

    maskdB = zeros(1,length(mask_sample_freqs_kHz));
    for m=1:length(masker_freqs_kHz)
      %maskdB = max(maskdB, schroeder(masker_freqs_kHz(m), mask_sample_freqs_kHz) + masker_amps_dB(m)); 
      maskdB = max(maskdB, parabolic_resonator(masker_freqs_kHz(m), mask_sample_freqs_kHz) + masker_amps_dB(m)); 
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


% Alternative mask function that uses parabolas for fast computetion.

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

function maskdB_ = decimate_frame_rate(model, decimate, f, frames, mask_sample_freqs_kHz)
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

    printf("f: %d left_f: %d right_f: %d left_fraction: %f right_fraction: %f \n", f, left_f, right_f, left_fraction, right_fraction)

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
    maskdB_left = interp1(left_sample_freqs_kHz, left_AmdB, sample_freqs_kHz);
    maskdB_right = interp1(right_sample_freqs_kHz, right_AmdB, sample_freqs_kHz);

    maskdB_ = left_fraction*maskdB_left + right_fraction*maskdB_right;
endfunction


% plot some masking curves, used for working on masking filter changes

function plot_masking
  Fs = 8000;

  figure(1)
  mask_sample_freqs_kHz = 0.1:0.025:(Fs/1000)/2;
  #{
  maskdB_s0 = schroeder(0.5, mask_sample_freqs_kHz, 0);
  plot(mask_sample_freqs_kHz, maskdB_s0,';schroeder 0;');
  maskdB_s1 = schroeder(0.5, mask_sample_freqs_kHz, 1);
  plot(mask_sample_freqs_kHz, maskdB_s1,'g;schroeder 1;');
  #}
  maskdB_res = parabolic_resonator(0.5, mask_sample_freqs_kHz);
  plot(mask_sample_freqs_kHz, maskdB_res,'r;resonator;');
  hold on;

  for f=0.5:0.5:3
    #{
    maskdB_s0 = schroeder(f, mask_sample_freqs_kHz, 0);
    plot(mask_sample_freqs_kHz, maskdB_s0);
    maskdB_s1 = schroeder(f, mask_sample_freqs_kHz, 1);
    plot(mask_sample_freqs_kHz, maskdB_s1,'g');
    #}
    maskdB_res = parabolic_resonator(f, mask_sample_freqs_kHz);
    plot(mask_sample_freqs_kHz, maskdB_res,'r;resonator;');
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


% determine 4 voicing bits based on 2 decimated voicing bits

function [v] = est_voicing_bits(v1, v5)
  if v1 == v5
    v(1:4) = v1;
  else
    v(1:2) = v1;
    v(3:4) = v5;
  end
endfunction


function [AmdB_ residual] = piecewise_model(AmdB, Wo) 
    L = length(AmdB);
    l1000 = floor(L/4);     
    AmdB_ = ones(1,L);
    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;

    % fit a resonator to max of first 300 - 1000 Hz

    fmin = 0.150;
    lmin = floor(L*fmin/4);
    [mx mx_ind] = max(AmdB(lmin+1:l1000));
    mx_ind += lmin;
    AmdB_ = parabolic_resonator(mx_ind*Wo*4/pi, mask_sample_freqs_kHz) + mx;  

    % fit a 2nd resonator, must be above 1000Hz

    fr1 = mx_ind*Wo*4/pi;
    fmin = 1;
    lmin = round(L*fmin/4);

    [mx mx_ind] = max(AmdB(lmin+1:L));
    mx_ind += lmin;
    AmdB_ = max(AmdB_, parabolic_resonator(mx_ind*Wo*4/pi, mask_sample_freqs_kHz) + mx);  
    fr2 = mx_ind*Wo*4/pi;

    % fit a third resonator, must be +/- 300 Hz after 2nd resonator

    residual = AmdB;
    residual(1:lmin) = 0;

    fr2 = mx_ind*Wo*4/pi;
    fmin = fr2 - 0.300;
    fmax = fr2 + 0.300;
    
    lmin = max(1, round(L*fmin/4));
    lmax = min(L, round(L*fmax/4));
    residual(lmin:lmax) = 0;

    [mx mx_ind] = max(residual);
    AmdB_ = max(AmdB_, parabolic_resonator(mx_ind*Wo*4/pi, mask_sample_freqs_kHz) + mx);  
    fr3 = mx_ind*Wo*4/pi;
   
    %figure(3);
    %plot(mask_sample_freqs_kHz, residual);

    printf("\nfr1: %f fr2: %f fr3: %f\n", fr1, fr2, fr3);
endfunction
