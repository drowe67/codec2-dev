% newamp1_batch.m
%
% Copyright David Rowe 2016
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Octave script to batch process model parameters using the new
% amplitude model.  Used for generating samples we can listen to.
%
% Usage:
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
%   $ cd ~/codec2-dev/octave
%   octave:14> newamp1_batch("../build_linux/src/hts1a")
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --amread hts1a_am.out -o - | play -t raw -r 8000 -s -2 -
% Or with a little more processing:
%   codec2-dev/build_linux/src$ ./c2sim ../../raw/hts2a.raw --amread hts2a_am.out --awread hts2a_aw.out --phase0 --postfilter --Woread hts2a_Wo.out -o - | play -q -t raw -r 8000 -s -2 -

% process a whole file and write results
% TODO: 
%   [ ] refactor decimate-in-time to avoid extra to/from model conversions
%   [ ] switches to turn on/off quantisation
%   [ ] rename mask_sample_freqs, do we need "mask" any more

#{
    [ ] refactor
#}

% In general, this function processes a bunch of amplitudes, we then
% use c2sim to hear the results

function surface = newamp1_batch(samname, optional_Am_out_name, optional_Aw_out_name)
  newamp;
  more off;

  max_amp = 80;
  postfilter = 0;

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);

  % Choose experiment to run test here -----------------------

  %model_ = experiment_filter(model);
  %model_ = experiment_filter_dec_filter(model);
  [model_ surface] = experiment_mel_freq(model, 0);
  %[model_ surface] = experiment_mel_diff_freq(model, 0);
  %model_ = experiment_dec_linear(model_);
  %[model_ rate_K_surface] = experiment_closed_loop_mean(model);

  % ----------------------------------------------------

  if nargin == 2
    Am_out_name = optional_Am_out_name;
  else
    Am_out_name = sprintf("%s_am.out", samname);
  end

  fam  = fopen(Am_out_name,"wb"); 

  Wo_out_name = sprintf("%s_Wo.out", samname);
  fWo  = fopen(Wo_out_name,"wb"); 
  
  for f=1:frames
    printf("%d ", f);   
    Wo = model_(f,1);
    L = min([model_(f,2) max_amp-1]);
    Am = model_(f,3:(L+2));

    Am_ = zeros(1,max_amp);
    Am_(2:L) = Am(1:L-1);

    % post filter, boosts higher amplitudes more than lower, improving
    % shape of formants and reducing muffling.  Note energy normalisation

    if postfilter
      e1 = sum(Am_(2:L).^2);
      Am_(2:L) = Am_(2:L) .^ 1.5;         
      e2 = sum(Am_(2:L).^2);
      Am_(2:L) *= sqrt(e1/e2);
    end

    fwrite(fam, Am_, "float32");
    fwrite(fWo, Wo, "float32");
  end
 
  fclose(fam);
  fclose(fWo);
  printf("\n")

endfunction
 
% Non linear sampling of frequency axis, reducing the "rate" is a
% first step before VQ

function mel = ftomel(fHz)
  mel = floor(2595*log10(1+fHz/700)+0.5);
endfunction


function [rate_K_surface rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model, K) 
  [frames nc] = size(model);
  mel_start = ftomel(200); mel_end = ftomel(3700); 
  step = (mel_end-mel_start)/(K-1);
  mel = mel_start:step:mel_end;
  rate_K_sample_freqs_Hz = 700*((10 .^ (mel/2595)) - 1);
  rate_K_sample_freqs_kHz = rate_K_sample_freqs_Hz/1000;

  rate_K_surface = resample_const_rate_f(model, rate_K_sample_freqs_kHz);
endfunction


function [model_ rate_K_surface] = experiment_mel_freq(model, vq_en=0)
  [frames nc] = size(model);
  K = 20; 
  [rate_K_surface  rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model, K);
  
  figure(1); clf; mesh(rate_K_surface);

  if vq_en
    melvq;
    load surface_vq; m=5;
   
    for f=1:frames
      mean_f(f) = mean(rate_K_surface(f,:));
      rate_K_surface_no_mean(f,:) = rate_K_surface(f,:) - mean_f(f);
    end
    
    [res rate_K_surface_ ind] = mbest(surface_vq, rate_K_surface_no_mean, m);

    % pf, needs some energy equalisation, does gd things for hts1a
    rate_K_surface_ *= 1.2;

    for f=1:frames
      rate_K_surface_(f,:) += mean_f(f);
    end

    rate_K_surface = rate_K_surface_;
     
    #{
    Nf = 4; Nf2 = 6;
    [b a]= cheby1(4, 1, 0.20);
    for k=1:K
      rate_K_surface_(:,k) = filter(b, a,  rate_K_surface_(:,k));
    end
    rate_K_surface_ = [rate_K_surface_(Nf2:frames,:); zeros(Nf2, K)];
    #}
  end

  model_ = resample_rate_L(model, rate_K_surface, rate_K_sample_freqs_kHz);

  %figure(2); clf; mesh(model_);

  
  for f=1:frames
    rate_K_surface(f,:) -= mean(rate_K_surface(f,:));
  end
endfunction


% mel spaced sampling, differential in time VQ.  Curiously, couldn't
% get very good results out of this, I suspect a bug

function [model_ rate_K_surface_diff] = experiment_mel_diff_freq(model, vq_en=0)
  [frames nc] = size(model);
  K = 20; 
  [rate_K_surface rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model, K);
  
  if vq_en
    melvq;
    load surface_diff_vq; m=5;
  end

  for f=1:frames
    mean_f(f,:) = mean(rate_K_surface(f,:));
    rate_K_surface_no_mean(f,:) = rate_K_surface(f,:) - mean_f(f,:);
  end

  rate_K_surface_no_mean_ = zeros(frames, K);
  rate_K_surface_no_mean_diff = zeros(frames, K);
  rate_K_surface_(1,:) = rate_K_surface_diff(1,:) = zeros(1, K);

  for f=2:frames
    rate_K_surface_diff(f,:) = rate_K_surface_no_mean(f,:) - 0.8*rate_K_surface_no_mean_(f-1,:);
    if vq_en
      [res arate_K_surface_diff_ ind] = mbest(surface_diff_vq, rate_K_surface_diff(f,:), m);
      rate_K_surface_diff_(f,:) = arate_K_surface_diff_;
    else
      rate_K_surface_diff_(f,:) = rate_K_surface_diff(f,:);
    end
    rate_K_surface_no_mean_(f,:) = 0.8*rate_K_surface_no_mean_(f-1,:) + rate_K_surface_diff_(f,:);
  end
  
  for f=1:frames
    rate_K_surface_(f,:) = rate_K_surface_no_mean_(f,:) + mean_f(f,:);
  end

  model_ = resample_rate_L(model, rate_K_surface_, rate_K_sample_freqs_kHz);

endfunction


% try vq with open and closed loop mean removal, turns out they give
% identical results lol

function [model_ rate_K_surface] = experiment_closed_loop_mean(model)
  [frames nc] = size(model);
  K = 15; 
  mel_start = ftomel(300); mel_end = ftomel(3000); 
  step = (mel_end-mel_start)/(K-1);
  mel = mel_start:step:mel_end;
  rate_K_sample_freqs_Hz = 700*((10 .^ (mel/2595)) - 1);
  rate_K_sample_freqs_kHz = rate_K_sample_freqs_Hz/1000;

  rate_K_surface = resample_const_rate_f(model, rate_K_sample_freqs_kHz);
  
  load surface_vq; m=1;
   
  for f=1:frames
    amean = mean(rate_K_surface(f,:));
    rate_K_target_no_mean = rate_K_surface(f,:) - amean;
    mse_open_loop(f) = search_vq2(surface_vq(:,:,1), rate_K_target_no_mean, m, 0);
    mse_closed_loop(f) = search_vq2(surface_vq(:,:,1), rate_K_surface(f,:), m, 1);
  end

  printf("rms open loop..: %f\nrms closed loop: %f\n", sqrt(mean(mse_open_loop)), sqrt(mean(mse_closed_loop)));

  % just return model_ as we have to so nothing breaks, it's not actually useful

  model_ = resample_rate_L(model, rate_K_surface, rate_K_sample_freqs_kHz);
    
endfunction


% conventional decimation in time without any filtering, then linear
% interpolation.  Linear interpolation is a two-tap (weak) form of fir
% filtering that may have problems with signals with high freq
% components, for example after quantisation noise is added.  Need to
% look into this some more.
 
function model_ = experiment_dec_linear(model)
  newamp;
  max_amp = 80;

  [frames nc] = size(model);
  model_ = zeros(frames, max_amp+3);
  decimate = 4;
  for f=1:frames    
    AmdB_ = decimate_frame_rate(model, decimate, f, frames);
    L = length(AmdB_);
    model_(f,1) = model(f,1); model_(f,2) = L; model_(f,3:(L+2)) = 10 .^ (AmdB_(1:L)/20);
  end
endfunction

#{ 
  Filtering time axis or surface, as a first step before decimation.
  So given surface, lets look at spectral content and see if we can
  reduce it while maintaining speech quality.  First step is to dft
  across time and plot.

  This just has one filtering step, which may help quantisation.  In practice
  we may need filtering before decimation and at the inerpolation stage
#}

function model_ =  experiment_filter(model)
  [frames nc] = size(model);
  K = 40; rate_K_sample_freqs_kHz = (1:K)*4/K;
  [rate_K_surface rate_K_sample_freqs_kHz] = resample_const_rate_f(model, rate_K_sample_freqs_kHz);
  
  Nf = 4; Nf2 = 6;
  [b a]= cheby1(4, 1, 0.20);
  %Nf = 20; Nf2 = 10;
  %b = fir1(Nf, 0.25); a = [1 zeros(1, Nf)];
  %Nf = 1; Nf2 = 1;
  %b = 1; a = [1 zeros(1, Nf)];

  %Nf = 2; Nf2 = 1;
  %beta = 0.99; w = pi/4;
  %b = [1 -2*beta*cos(w) beta*beta]; a = [1 zeros(1, Nf)];
  %Nf = 10; Nf2 = 10;
  %b = fir2(10, [0 0.2 0.3 1], [1 1 0.1 0.1]); a = [1 zeros(1, Nf)];

  %Nf = 1; Nf2 = 1;
  dft_surface = zeros(frames,K);
  rate_K_surface_filt = zeros(frames,K);
  dft_surface_filt = zeros(frames,K);
  for k=1:K
    dft_surface(:,k) = fft(rate_K_surface(:,k).*hanning(frames));
    %rate_K_surface_filt(:,k) = filter(b, a, rate_K_surface(:,k));
    rate_K_surface_filt(:,k) = filter(b, a,  rate_K_surface(:,k));
    dft_surface_filt(:,k) = fft(rate_K_surface_filt(:,k).*hanning(frames));
  end
  figure(1); clf;
  mesh(abs(dft_surface))
  figure(2); clf;
  mesh(abs(dft_surface_filt))

  Fs = 100; Ts = 1/Fs;
  figure(3);
  subplot(211);
  h = freqz(b,a,Fs/2);
  plot(1:Fs/2, 20*log10(abs(h)))
  axis([1 Fs/2 -40 0])
  ylabel('Gain (dB)')
  grid;
  subplot(212)
  [g w] = grpdelay(b,a);
  plot(w*Fs/pi, 1000*g*Ts)
  %axis([1 Fs/2 0 0.5])
  xlabel('Frequency (Hz)');
  ylabel('Delay (ms)')
  grid

  figure(4); clf; stem(filter(b,a,[1 zeros(1,20)]))
  
  % adjust for time offset due to filtering

  rate_K_surface_filt = [rate_K_surface_filt(Nf2:frames,:); zeros(Nf2, K)];

  % back down to rate L

  model_ = resample_rate_L(model, rate_K_surface_filt, rate_K_sample_freqs_kHz);
endfunction


% filter, decimate, zero insert, filter, simulates decimation and re-interpolation, as
% an alternative to simple linear interpolation.

function model_ =  experiment_filter_dec_filter(model)
  [frames nc] = size(model);

  % rate K surface

  K = 40; rate_K_sample_freqs_kHz = (1:K)*4/K;
  [rate_K_surface rate_K_sample_freqs_kHz] = resample_const_rate_f(model, rate_K_sample_freqs_kHz);
  
  % filter, not we run across each of the K bins, treating them as a 1-D sequence

  Nf = 4; filter_delay = 12;
  [b a]= cheby1(4, 1, 0.20);
  %Nf = 10; Nf2 = Nf/2; filter_delay = Nf2*2;
  %b = fir1(Nf, 0.25); a = [1 zeros(1, Nf)];
  %Nf = 10; Nf2 = 2; filter_delay = 2*Nf2;
  %b = fir2(10, [0 0.2 0.3 1], [1 1 0.1 0.1]); a = [1 zeros(1, Nf)];

  rate_K_surface_filt = zeros(frames,K);
  for k=1:K
    rate_K_surface_filt(:,k) = filter(b, a, rate_K_surface(:,k));
  end

  % decimate from 100 to 25Hz, and zero pad, which we simulate by
  % setting all K samples in 3 out of 4 time-samples to zero

  M = 4;

  for f=1:frames
    if mod(f,M)
      rate_K_surface_filt(f,:) = zeros(1,K);
    end
  end

  % filter to reconstruct 100 Hz frame rate

  rate_K_surface_filt_recon = zeros(frames,K);
  for k=1:K
    rate_K_surface_filt_recon(:,k) = filter(b, a, M*rate_K_surface_filt(:,k));
  end

  figure(1); clf;
  mesh(rate_K_surface_filt)
  figure(2); clf;
  mesh(rate_K_surface_filt_recon)

  % adjust for time offset due to 2 lots of filtering

  rate_K_surface_filt_recon = [rate_K_surface_filt_recon(filter_delay:frames,:); zeros(filter_delay, K)];

  % back down to rate L

  model_ = resample_rate_L(model, rate_K_surface_filt_recon, rate_K_sample_freqs_kHz);
endfunction


#{
todo: get this working again

function model_ = experiment_piecewise(model)
  % encoder loop ------------------------------------------------------

  fvec_log = []; amps_log = [];

  for f=1:frames
    printf("%d ", f);   
    Wo = model(f,1);
    L = min([model(f,2) max_amp-1]);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);
    e(f) = sum(AmdB)/L;

    % fit model

    [AmdB_ res fvec fvec_ amps] = piecewise_model(AmdB, Wo);
    fvec_log = [fvec_log; fvec];
    amps_log = [amps_log; amps];

    model_(f,1) = Wo; model_(f,2) = L; model_(f,3:(L+2)) = 10 .^ (AmdB(1:L)/20);
  end

  for f=1:frames
    if 0
    if (f > 1) && (e(f) > (e(f-1)+3))
        decimate = 2;
      else
        decimate = 4;
    end
    end
    printf("%d ", decimate);   
    
    AmdB_ = decimate_frame_rate(model_, decimate, f, frames);
    L = length(AmdB_);

    Am_ = zeros(1,max_amp);
    Am_(2:L) = 10 .^ (AmdB_(1:L-1)/20);  % C array doesnt use A[0]
    fwrite(fam, Am_, "float32");
    fwrite(fWo, Wo, "float32");
  end
endfunction

#}



% Resample Am from time-varying rate L=floor(pi/Wo) to fixed rate K.  This can be viewed
% as a 3D surface with time, freq, and ampitude axis.

function [rate_K_surface rate_K_sample_freqs_kHz] = resample_const_rate_f(model, K=50)

  % convert rate L=pi/Wo amplitude samples to fixed rate K

  max_amp = 80;
  [frames col] = size(model);
  rate_K_sample_freqs_kHz = (1:K)*4/K;
  rate_K_surface = zeros(frames, K);

  for f=1:frames
    Wo = model(f,1);
    L = min([model(f,2) max_amp-1]);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);
    rate_L_sample_freqs_kHz = (1:L)*Wo*4/pi;
    
    rate_K_surface(f,:) = interp1(rate_L_sample_freqs_kHz, AmdB, rate_K_sample_freqs_kHz, "spline", "extrap");
  end
endfunction


function [rate_K_surface rate_K_sample_freqs_kHz] = resample_const_rate_f(model, rate_K_sample_freqs_kHz)

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
    rate_L_sample_freqs_kHz = (1:L)*Wo*4/pi;
    
    rate_K_surface(f,:) = interp1(rate_L_sample_freqs_kHz, AmdB, rate_K_sample_freqs_kHz, "spline", 0);
  end
endfunction


% Take a rate K surface and convert back to time varying rate L

function model_ = resample_rate_L(model, rate_K_surface, rate_K_sample_freqs_kHz)
  max_amp = 80;
  [frames col] = size(model);

  model_ = zeros(frames, max_amp+3);
  for f=1:frames-1
    Wo = model(f,1);
    L = min(pi/Wo, max_amp-1);
    rate_L_sample_freqs_kHz = (1:L)*Wo*4/pi;
    
    % back down to rate L

    AmdB_ = interp1(rate_K_sample_freqs_kHz, rate_K_surface(f,:), rate_L_sample_freqs_kHz, "spline", 0);

    model_(f,1) = Wo; model_(f,2) = L; model_(f,3:(L+2)) = 10 .^ (AmdB_(1:L)/20);
   end
endfunction


% early test, devised to test rate K<->L changes along frequency axis

function model_ = resample_half_frame_offset(model, rate_K_surface, rate_K_sample_freqs_kHz)
  max_amp = 80;
  [frames col] = size(model);

  % Check distortion by shifting half a frame and returning model for
  % synthesis Hmm how to handle Wo?  Well lets assume it's OK and
  % average it.  If it sounds OK then it must be OK as well and Am
  % interpolation.  If not then will do some more thinking.

  model_ = zeros(frames, max_amp+3);
  for f=1:frames-1
    Wo = (model(f,1) + model(f+1,1))/2;
    L = min(pi/Wo, max_amp-1);
    rate_L_sample_freqs_kHz = (1:L)*Wo*4/pi;
    
    % interpolate at half frame offset

    half = 0.5*rate_K_surface(f,:) + 0.5*rate_K_surface(f+1,:);
    
    % back down to rate L

    AmdB_ = interp1(rate_K_sample_freqs_kHz, half, rate_L_sample_freqs_kHz, "spline", "extrap");

    % todo a way to save all model params ... Am, Wo, L, (v?)  Start with current facility then
    % make it better

    model_(f,1) = Wo; model_(f,2) = L; model_(f,3:(L+2)) = 10 .^ (AmdB_(1:L)/20);
   end
endfunction


% vq search with optional closed loop mean estimation, turns out this gives
% identical results to extracting the mean externally

function [mse_list index_list] = search_vq2(vq, target, m, closed_loop_dc = 0)

  [Nvec order] = size(vq);

  mse = zeros(1, Nvec);

  % find mse for each vector

  for i=1:Nvec
     if closed_loop_dc
       sum(target - vq(i,:))
       g = sum(target - vq(i,:))/order;
       mse(i) = sum((target - vq(i,:) - g) .^2);
     else
       mse(i) = sum((target - vq(i,:)) .^2);
    end
  end

  % sort and keep top m matches

  [mse_list index_list ] = sort(mse);

  mse_list = mse_list(1:m);
  index_list = index_list(1:m);

endfunction
