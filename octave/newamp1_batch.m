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

function [fvec_log amps_log] = newamp1_batch(samname, optional_Am_out_name, optional_Aw_out_name)
  newamp;
  more off;

  max_amp = 80;

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);

  model_ = experiment_filter(model);

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
    fwrite(fam, Am_, "float32");
    fwrite(fWo, Wo, "float32");
  end
 
  fclose(fam);
  fclose(fWo);
  printf("\n")

endfunction
 
 
#{ 
  Filtering time axis or surface, as a first step before decimation.
  So given surface, lets look at spectral content and see if we can
  reduce it while maintaining speech quality.  First step is to dft
  across time and plot.
#}

function model_ =  experiment_filter(model)
  [frames nc] = size(model);
  K = 40;
  [rate_K_surface rate_K_sample_freqs_kHz] = resample_const_rate_f(model, K);
  
  Nf = 4; Nf2 = 5;
  [b a]= cheby1(4, 1, 0.20);
  dft_surface = zeros(frames,K);
  rate_K_surface_filt = zeros(frames,K);
  dft_surface_filt = zeros(frames,K);
  for k=1:K
    dft_surface(:,k) = fft(rate_K_surface(:,k).*hanning(frames));
    rate_K_surface_filt(:,k) = filter(b, a, rate_K_surface(:,k));
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
  axis([1 Fs/2 -80 0])
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
% as a 3D surface with tim, freq, nd ampitude axis.


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


% Take a rate K surface and convert back to time varying rate L

function model_ = resample_rate_L(model, rate_K_surface, rate_K_sample_freqs_kHz, K=50)
  max_amp = 80;
  [frames col] = size(model);

  model_ = zeros(frames, max_amp+3);
  for f=1:frames-1
    Wo = (model(f,1) + model(f+1,1))/2;
    L = min(pi/Wo, max_amp-1);
    rate_L_sample_freqs_kHz = (1:L)*Wo*4/pi;
    
    % back down to rate L

    AmdB_ = interp1(rate_K_sample_freqs_kHz, rate_K_surface(f,:), rate_L_sample_freqs_kHz, "spline", "extrap");

    % todo a way to save all model params ... Am, Wo, L, (v?)  Start with current facility then
    % make it better

    model_(f,1) = Wo; model_(f,2) = L; model_(f,3:(L+2)) = 10 .^ (AmdB_(1:L)/20);
   end
endfunction


function model_ = resample_half_frame_offset(model, rate_K_surface, rate_K_sample_freqs_kHz, K=50)
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
