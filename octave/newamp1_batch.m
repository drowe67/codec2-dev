% newamp1_batch.m
%
% Copyright David Rowe 2016
% This program is distributed under the terms of the GNU General Public License 
% Version 2

#{

  Octave script to batch process model parameters using the new
  amplitude model, version 1.  Outputs another set of model parameters
  that can be fed to c2sim for listening tests.  The companion
  newamp1_fbf.m script is used to visualise the processing frame by frame
 
  c2sim -> dump files -> newamp1_batch.m -> output model params -> c2sim -> play
 
  The newamp1_xxx scripts have evolved to (i) resample {Am} using a
  mel frequency axis, (ii) 2 stage VQ the mean removed vector.  Seems to work
  OK at 700 bit/s, comparable to 1300.

  Usage:

    build codec2 with -DDUMP - see codec2-dev/README, then:

    ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
    $ cd ~/codec2-dev/octave
    octave:14> newamp1_batch("../build_linux/src/hts1a")
    ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --amread hts1a_am.out -o - | play -t raw -r 8000 -s -2 -
 
  Or with a little more processing, first dump energy and voicing, the
  import Wo, voicing, phase spectra which simulates all the decoder
  DSP.  We switch on lpc 10 just to dump voicing.

    $ ./c2sim ../../raw/vk5qi.raw --phase0 --postfilter --dump vk5qi --lpc 10 --dump_pitch_e vk5qi_pitche.txt
    octave:14> newamp1_batch("../build_linux/src/vk5qi", "../build_linux/src/vk5qi.out");
    $ ./c2sim ../../raw/vk5qi.raw --phase0 --postfilter --amread vk5qi_am.out --hmread vk5qi_hm.out --Woread vk5qi_Wo.out --hand_voicing vk5qi_v.txt -o - | play -q -t raw -r 8000 -s -2 -

#}


% In general, this function processes a bunch of amplitudes, we then
% use c2sim to hear the results.  Bunch of different experiments below

function [surface mean_f] = newamp1_batch(input_prefix, varargin)
  newamp;
  more off;

  max_amp = 160;
  mean_f = [];

  % defaults

  synth_phase = output = 1;
  output_prefix = input_prefix;
  vq_type = "";
  vq_filename = "";
  fit_order = 0;
  mode = "const";

  % parse variable argument list

  if (length (varargin) > 0)

    % check for the "output_prefix" option

    ind = arg_exists(varargin, "output_prefix");
    if ind
      output_prefix =  varargin{ind+1};
    end
    ind = arg_exists(varargin, "mode");
    if ind
      mode =  varargin{ind+1};
    end
    ind = arg_exists(varargin, "nomean");
    if ind
      fit_order = 0;
    end
    ind = arg_exists(varargin, "noslope");
    if ind
      fit_order = 1;
    end
    ind = arg_exists(varargin, "vqh");
    if ind
      vq_type = varargin{ind};
      vq_filename = varargin{ind+1};
      printf("vq_type: %s vq_filename: %s\n", vq_type, vq_filename);
    end
    ind = arg_exists(varargin, "no_output");
    if ind
      output = 0;
      synth_phase = 0;
    end
  end

  printf("output_prefix: %s\nfit_order %d output: %d\n", output_prefix, fit_order, output);
  
  model_name = strcat(input_prefix,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);

  voicing_name = strcat(input_prefix,"_pitche.txt");
  voicing = zeros(1,frames);
  
  if exist(voicing_name, "file") == 2
    pitche = load(voicing_name);
    voicing = pitche(:, 3);
  end

  % Choose experiment to run test here -----------------------

  if strcmp(mode, 'dct2')
    [model_ surface] = experiment_rate_K_dct2(model, 0, 1, voicing);
  end
  if strcmp(mode, 'mel')
    [model_ surface] = experiment_mel_freq(model, 0, 1, voicing);
  end
  if strcmp(mode, 'const')
    [model_ surface] = experiment_const_freq(model, fit_order, vq_type, vq_filename);
  end
  if strcmp(mode, 'piecewise')
    model_ = experiment_piecewise(model);
  end
  if strcmp(mode, 'dct')
    model_ = experiment_dct(model);
  end

  % ----------------------------------------------------

  if output
    Am_out_name = sprintf("%s_am.out", output_prefix);
    fam  = fopen(Am_out_name,"wb"); 
 
    Wo_out_name = sprintf("%s_Wo.out", output_prefix);
    fWo  = fopen(Wo_out_name,"wb"); 
  
    if synth_phase
     Hm_out_name = sprintf("%s_hm.out", output_prefix);
     fhm = fopen(Hm_out_name,"wb"); 
    end

    for f=1:frames
      %printf("%d ", f);   
      Wo = model_(f,1); L = min([model_(f,2) max_amp-1]); Am = model_(f,3:(L+2));
      if Wo*L > pi
        printf("Problem: %d  Wo*L > pi\n", f);   
      end

      Am_ = zeros(1,max_amp); Am_(2:L) = Am(1:L-1); fwrite(fam, Am_, "float32");
      fwrite(fWo, Wo, "float32");

      if synth_phase

        % synthesis phase spectra from magnitiude spectra using minimum phase techniques

        fft_enc = 256;
        phase = determine_phase(model_, f, fft_enc);
        assert(length(phase) == fft_enc);

        % sample phase at centre of each harmonic, not 1st entry Hm[1] in octave Hm[0] in C
        % is not used

        Hm = zeros(1, 2*max_amp);
        for m=1:L
          b = round(m*Wo*fft_enc/(2*pi));
          Hm(2*m) = cos(phase(b));
          Hm(2*m+1) = -sin(phase(b));
        end
        fwrite(fhm, Hm, "float32");    
      end
    end
 
   fclose(fam);
   fclose(fWo);
   if synth_phase
     fclose(fhm);
   end

   % save voicing file
  
   if exist("voicing_", "var")
     v_out_name = sprintf("%s_v.txt", output_prefix);
     fv  = fopen(v_out_name,"wt"); 
     for f=1:length(voicing_)
       fprintf(fv,"%d\n", voicing_(f));
     end
     fclose(fv);
    end
  end

  if fit_order == 0;
    for f=1:frames
      surface(f,:) -= mean(surface(f,:));
    end
  end

  if fit_order == 1
    surface = slope_and_mean_removal(surface);
  end

  printf("\n")

endfunction
 

function surface = slope_and_mean_removal(surface)
  [frames K] = size(surface);
  for f=1:frames
    v = surface(f,:);
    [m b] = linreg(1:K, v, K);
    v -= m*(1:K) + b;
    surface(f,:) = v;
  end
endfunction


function ind = arg_exists(v, str) 
   ind = 0;
   for i=1:length(v)
      if strcmp(v{i}, str)
        ind = i;
      end     
    end
endfunction


% Basic unquantised rate K linear sampling then back to rate L.  Used for generating 
% training vectors and testing vector quntisers.

function [model_ rate_K_surface] = experiment_const_freq(model, fit_order, vq_type, vq_filename)
  melvq;

  [frames nc] = size(model);
  Fs = 8000;
  fg = 1;
  quant_en = 0;

  rate_K_sample_freqs_kHz = [0.1:0.1:4];
  K = length(rate_K_sample_freqs_kHz);

  if strcmp(vq_type, "vql")
    quant_en = 1;
    x = load(vq_filename); vq = x.vq;
    [vq_rows vq_cols] = size(vq); vq_st = 1; vq_en = vq_st + vq_cols - 1;
  end

  if strcmp(vq_type, "vqh")
    quant_en = 1;
    x = load(vq_filename); vq = x.vq;
    [vq_rows vq_cols] = size(vq); vq_st = 11; vq_en = K;
  end

  energy = zeros(1,frames);
  for f=1:frames
    L = model(f,2);
    energy(f) = 10*log10(sum( model(f,3:(L+2)) .^ 2 ));
  end

  rate_K_surface = resample_const_rate_f(model, rate_K_sample_freqs_kHz, Fs);

  rate_K_surface_fit = zeros(frames,K);
  b = slope = zeros(1,frames);

  if fit_order == 0
    for f=1:frames
      b(f) = mean(rate_K_surface(f,:));
      rate_K_surface_fit(f,:) = rate_K_surface(f,:) - b(f);
    end
  end

  if fit_order == 1
    for f=1:frames
      [aslope ab] = linreg(1:K, rate_K_surface(f,:), K);
      rate_K_surface_fit(f,:) = rate_K_surface(f,:) - aslope*(1:K) - ab;
      slope(f) = aslope; b(f) = ab;
    end
  end

  % optional vector quantise
  
  if quant_en
 
    rate_K_surface_fit_ = zeros(frames, K);
    res = zeros(frames,vq_cols); ind = [];

    for f=1:frames
      target = rate_K_surface_fit(f, vq_st:vq_en);

      [diff_weighted weights error g mn_ind] = search_vq_weighted(target, vq);
      if (f>=73) && (f<=75)
        printf("f: %d mn_ind: %d g: %3.2f sd: %3.2f\n", f, mn_ind, g(mn_ind), error(mn_ind));
      end
      rate_K_surfacurface_fit_(f, vq_st:vq_en) = vq(mn_ind,:) + g(mn_ind);

      %res(f,vq_st:vq_en) = rate_K_surface_no_mean(f,vq_st:vq_en) - rate_K_surface_no_mean_(f,vq_st:vq_en);
      res(f,vq_st:vq_en) = diff_weighted(mn_ind,:);
      ind = [ind mn_ind];
    end
    
    figure(fg++); clf; mesh(res);
    sd_per_frame = std(res(:,vq_st:vq_en)');
    figure(fg++); subplot(211); plot(energy); subplot(212); plot(sd_per_frame);
    figure(fg++); subplot(211); hist(sd_per_frame); subplot(212); hist(ind,100);
    printf("VQ rms SD: %3.2f\n", mean(sd_per_frame));
  else
    rate_K_surface_fit_ = rate_K_surface_fit;
  end

  for f=1:frames
    rate_K_surface_(f,:) = rate_K_surface_fit_(f,:) + slope(f)*(1:K) + b(f);
  end

  [model_ AmdB_] = resample_rate_L(model, rate_K_surface_, rate_K_sample_freqs_kHz, Fs);

  % Measure distortion between AmdB and AmdB_, ths includes distortion
  % in the rate K <-> L transition.  Can optionally plot distorted
  % frames

  plot_sd_thresh = 10;
  sd = zeros(1,frames);
  for f=1:frames
    Wo = model(f,1);
    L = model(f,2);
    AmdB = 20*log10(model(f,3:(L+2)));
    sd(f) = std(AmdB(1:L) - AmdB_(f,1:L));
    if (sd(f) > plot_sd_thresh) && (fg < 10)
      printf("fg: %d f: %d\n", fg, f);
      figure(fg++); clf; plot((1:L)*Wo*4/pi, AmdB(1:L),'b+-'); hold on; plot((1:L)*Wo*4/pi, AmdB_(f,1:L),'r+-');
      plot(rate_K_sample_freqs_kHz, rate_K_surface_(f,:), 'c+-'); hold off;
    end
  end
  printf("rate K resampling SD: %3.2f\n", mean(sd));
  figure(fg++); clf; subplot(211); plot(energy); subplot(212); plot(sd);
  figure(fg++); clf; hist(sd);
endfunction


function model_ = experiment_dct(model)

  [frames tmp] = size(model); max_amp = 160;

  for f=1:frames
    printf("%d ", f);   
    Wo = model(f,1);
    L = min([model(f,2) max_amp-1]);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);

    % fit model

    D = dct(AmdB);
    E = zeros(1,L);
    E(1:min(20,L)) = D(1:min(20,L));
    AmdB_ = idct(D);

    model_(f,1) = Wo; model_(f,2) = L; model_(f,3:(L+2)) = 10 .^ (AmdB_(1:L)/20);
  end

endfunction


% -----------------------------------------------------------------------------------------
% Linear decimator/interpolator that operates at rate K, includes VQ, post filter, and Wo/E
% quantisation.  Evolved from abys decimator below.  Simulates the entire encoder/decoder.

function [model_ voicing_ indexes] = experiment_rate_K_dec(model, voicing)
  max_amp = 80;
  [frames nc] = size(model);
  model_ = zeros(frames, max_amp+3);
  indexes = zeros(frames,4);

  M = 4;
 
  % create frames x K surface.  TODO make all of this operate frame by
  % frame, or at least M/2=4 frames rather than one big chunk

  K = 20; 
  [surface sample_freqs_kHz] = resample_const_rate_f_mel(model, K);
  target_surface = surface;

  figure(1);
  mesh(surface);

  % VQ rate K surface.  TODO: If we are decimating by M/2=4 we really
  % only need to do this every 4th frame.

  melvq;
  load train_120_vq; m=5;
       
  for f=1:frames
    mean_f(f) = mean(surface(f,:));
    surface_no_mean(f,:) = surface(f,:) - mean_f(f);
  end
  figure(2);
  mesh(surface_no_mean);

  [res surface_no_mean_ ind] = mbest(train_120_vq, surface_no_mean, m);
  indexes(:,1:2) = ind;

  for f=1:frames
    surface_no_mean_(f,:) = post_filter(surface_no_mean_(f,:), sample_freqs_kHz, 1.5);
  end
  figure(3);
  mesh(surface_no_mean_);
    
  surface_ = zeros(frames, K);
  energy_q = create_energy_q;
  for f=1:frames   
    [mean_f_ indx] = quantise(energy_q, mean_f(f));
    indexes(f,3) = indx - 1;
    %mean_f_ = mean_f(f);
    surface_(f,:) = surface_no_mean_(f,:) + mean_f_;
  end

  figure();
  mesh(surface_);

  % break into segments of M frames.  We have 3 samples in M frame
  % segment spaced M/2 apart and interpolate the rest.  This evolved
  % from AbyS scheme below but could be simplified to simple linear
  % interpolation, or using 3 or 4 points but shift of M/2=4 frames.
  
  interpolated_surface_ = zeros(frames, K);
  for f=1:M:frames-M
    left_vec = surface_(f,:);
    right_vec = surface_(f+M,:);
    sample_points = [f f+M];
    resample_points = f:f+M-1;
    for k=1:K
      interpolated_surface_(resample_points,k) = interp_linear(sample_points, [left_vec(k) right_vec(k)], resample_points);
    end    
  end

  % break into segments for purposes of Wo interpolation

  for f=1:M:frames
    % quantise Wo

    % UV/V flag is coded using a zero index for Wo, this means we need to
    % adjust Wo index slightly for the lowest Wo V frames

    if voicing(f)
      index = encode_log_Wo(model(f,1), 6);
      if index == 0
        index = 1;
      end
      indexes(f,4) = index;
      model_(f,1) = decode_log_Wo(indexes(f,4), 6);
    else
      indexes(f,4) = 0;
      model_(f,1) = 2*pi/100;
    end
  end      


  voicing_ = zeros(1, frames);
  for f=1:M:frames-M

    Wo1_ = model_(f,1);
    Wo2_ = model_(f+M,1);

    % uncomment to use unquantised values
    %Wo1_ = model(f,1);
    %Wo2_ = model(f+M,1);

    if !voicing(f) && !voicing(f+M)
       model_(f:f+M-1,1) = 2*pi/100;
    end

    if voicing(f) && !voicing(f+M)
       model_(f:f+M/2-1,1) = Wo1_;
       model_(f+M/2:f+M-1,1) = 2*pi/100;
       voicing_(f:f+M/2-1) = 1;
    end

    if !voicing(f) && voicing(f+M)
       model_(f:f+M/2-1,1) = 2*pi/100;
       model_(f+M/2:f+M-1,1) = Wo2_;
       voicing_(f+M/2:f+M-1) = 1;
    end

    if voicing(f) && voicing(f+M)
      Wo_samples = [Wo1_ Wo2_];
      model_(f:f+M-1,1) = interp1([f f+M], Wo_samples, f:f+M-1, "linear", 0);
      voicing_(f:f+M-1) = 1;
    end

    #{
    printf("f: %d f+M/2: %d Wo: %f %f (%f %%) v: %d %d \n", f, f+M/2, model(f,1), model(f+M/2,1), 100*abs(model(f,1) - model(f+M/2,1))/model(f,1), voicing(f), voicing(f+M/2));
    for i=f:f+M/2-1
      printf("  f: %d v: %d v_: %d Wo: %f Wo_: %f\n", i, voicing(i), voicing_(i), model(i,1),  model_(i,1));
    end
    #}
  end
  model_(frames-M:frames,1) = pi/100; % set end frames to something sensible

  % enable these to use original (non interpolated) voicing and Wo
  %voicing_ = voicing;
  %model_(:,1) = model(:,1);

  model_(:,2) = floor(pi ./ model_(:,1)); % calculate L for each interpolated Wo
  model_ = resample_rate_L(model_, interpolated_surface_, sample_freqs_kHz);

endfunction


% ---------------------------------------------------------------------------------------
% Stand alone decoder that takes indexes and creates model_, models
% decoder and an important step in proving everything works
 
function [model_ voicing_] = model_from_indexes(indexes)
  max_amp = 80;  K = 20;  M = 4;

  [frames nc] = size(indexes);
  model = model_ = zeros(frames, max_amp+3);
  sample_freqs_kHz = mel_sample_freqs_kHz(K);
  energy_q = 10 + 40/16*(0:15);

  melvq;
  load train_120_vq;

  % decode vector quantised surface

  surface_no_mean_ = zeros(frames,K);
  surface_ = zeros(frames, K);
  for f=1:M:frames
    surface_no_mean_(f,:) = train_120_vq(indexes(f,1),:,1) + train_120_vq(indexes(f,2),:,2);
    surface_no_mean_(f,:) = post_filter(surface_no_mean_(f,:), sample_freqs_kHz, 1.5);
    mean_f_ = energy_q(indexes(f,3)+1);
    surface_(f,:) = surface_no_mean_(f,:) + mean_f_;
  end
    
  % break into segments of M frames.  We have 2 samples spaced M apart
  % and interpolate the rest.
  
  interpolated_surface_ = zeros(frames, K);
  for f=1:M:frames-M
    left_vec = surface_(f,:);
    right_vec = surface_(f+M,:);
    sample_points = [f f+M];
    resample_points = f:f+M-1;
    for k=1:K
      interpolated_surface_(resample_points,k) = interp_linear(sample_points, [left_vec(k) right_vec(k)], resample_points);
    end    
  end

  % recover Wo and voicing

  voicing = zeros(1, frames);
  for f=1:M:frames
    if indexes(f,4) == 0
      voicing(f) = 0;
      model(f,1) = 2*pi/100;
    else
      voicing(f) = 1;
      model(f,1) = decode_log_Wo(indexes(f,4), 6);
    end
  end

  % break into M segments for purposes of Wo interpolation

  voicing_ = zeros(1, frames);
  for f=1:M:frames-M

    Wo1_ = model(f,1);
    Wo2_ = model(f+M,1);

    if !voicing(f) && !voicing(f+M)
       model_(f:f+M-1,1) = 2*pi/100;
    end

    if voicing(f) && !voicing(f+M)
       model_(f:f+M/2-1,1) = Wo1_;
       model_(f+M/2:f+M-1,1) = 2*pi/100;
       voicing_(f:f+M/2-1) = 1;
    end

    if !voicing(f) && voicing(f+M)
       model_(f:f+M/2-1,1) = 2*pi/100;
       model_(f+M/2:f+M-1,1) = Wo2_;
       voicing_(f+M/2:f+M-1) = 1;
    end

    if voicing(f) && voicing(f+M)
      Wo_samples = [Wo1_ Wo2_];
      model_(f:f+M-1,1) = interp_linear([f f+M], Wo_samples, f:f+M-1);
      voicing_(f:f+M-1) = 1;
    end

    #{
    printf("f: %d f+M/2: %d Wo: %f %f (%f %%) v: %d %d \n", f, f+M/2, model(f,1), model(f+M/2,1), 100*abs(model(f,1) - model(f+M/2,1))/model(f,1), voicing(f), voicing(f+M/2));
    for i=f:f+M/2-1
      printf("  f: %d v: %d v_: %d Wo: %f Wo_: %f\n", i, voicing(i), voicing_(i), model(i,1),  model_(i,1));
    end
    #}
  end
  model_(frames-M:frames,1) = pi/100; % set end frames to something sensible

  % enable these to use original (non interpolated) voicing and Wo
  %voicing_ = voicing;
  %model_(:,1) = model(:,1);

  model_(:,2) = floor(pi ./ model_(:,1)); % calculate L for each interpolated Wo
  model_ = resample_rate_L(model_, interpolated_surface_, sample_freqs_kHz);

endfunction


% ---------------------------------------------------------------------------------------
% Stand alone decoder that takes indexes and creates model_, just like 
% model_from_indexes above. This version is refactored to perform frame by frame
% processing, as a stepping stone to C.

function [model_ voicing_] = model_from_indexes_fbf(indexes)
  max_amp = 80;  K = 20;  M = 4;

  [frames nc] = size(indexes);
  model = model_ = zeros(frames, max_amp+3);
  sample_freqs_kHz = mel_sample_freqs_kHz(K);
  energy_q = 10 + 40/16*(0:15);

  melvq;
  load train_120_vq;

  surface_no_mean_ = zeros(frames,K);
  surface_ = zeros(frames, K);
  interpolated_surface_ = zeros(frames, K);
  voicing = zeros(1, frames);
  voicing_ = zeros(1, frames);

  for f=1:M:frames
    % decode vector quantised surface

    surface_no_mean_(f,:) = train_120_vq(indexes(f,1),:,1) + train_120_vq(indexes(f,2),:,2);
    surface_no_mean_(f,:) = post_filter(surface_no_mean_(f,:), sample_freqs_kHz, 1.5);
    mean_f_ = energy_q(indexes(f,3)+1);
    surface_(f,:) = surface_no_mean_(f,:) + mean_f_;

    % break into segments of M frames.  We have 2 samples spaced M apart
    % and interpolate the rest.

    if f > M
      left_vec = surface_(f-M,:);
      right_vec = surface_(f,:);
      sample_points = [f-M f];
      resample_points = f-M:f-1;
      for k=1:K
        interpolated_surface_(resample_points,k) = interp_linear(sample_points, [left_vec(k) right_vec(k)], resample_points);
      end    
    end

    % recover Wo and voicing

    if indexes(f,4) == 0
      voicing(f) = 0;
      model(f,1) = 2*pi/100;
    else
      voicing(f) = 1;
      model(f,1) = decode_log_Wo(indexes(f,4), 6);
    end

    if f > M
      Wo1 = model(f-M,1);
      Wo2 = model(f,1);

      [Wo_ avoicing_] = interp_Wo_v(Wo1, Wo2, voicing(f-M), voicing(f));
      model_(f-M:f-1,1) = Wo_;
      voicing_(f-M:f-1) = avoicing_;
      model_(f-M:f-1,2) = floor(pi ./ model_(f-M:f-1,1)); % calculate L for each interpolated Wo
    end

  end
  
  model_(frames-M:frames,1) = pi/100; % set end frames to something sensible
  model_(frames-M:frames,2) = floor(pi ./ model_(frames-M:frames,1));

  model_ = resample_rate_L(model_, interpolated_surface_, sample_freqs_kHz);

endfunction



% ---------------------------------------------------------------------------------------
% Various experiments tried during development

% rate K mel-resampling, high end correction, and DCT experiment workhorse

function [model_ rate_K_surface] = experiment_rate_K_dct2(model, vq_en=0, plots=1, voicing)
  [frames nc] = size(model);
  K = 20; Fs = 8000; correct_rate_K_en = 1;

  %quantisers = load("dct2quant.txt"); 
  %nlevels = load("dct2quantnlevels.txt");
  %map = load("dct2map.txt");

  for f=1:frames
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;
    [rate_K_vec rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model(f,:), K);
    if correct_rate_K_en
      [tmp_ AmdB_] = resample_rate_L(model(f,:), rate_K_vec, rate_K_sample_freqs_kHz);
      [rate_K_vec_corrected orig_error error nasty_error_log nasty_error_m_log] = correct_rate_K_vec(rate_K_vec, rate_K_sample_freqs_kHz, AmdB, AmdB_, K, Wo, L, Fs);
      rate_K_surface(f,:) = rate_K_vec_corrected;
    else
      rate_K_surface(f,:) = rate_K_vec;
    end
  end

  % break into 160ms blocks, 2D DCT, truncate, IDCT

  Tf = 0.01;   % frame period in seconds
  Nt = 16;     % number of 10ms frames blocks in time
  dec = 2;     % decimation factor
  dist_dB = 2; % use enough coefficients to get this distortion ond DCT coeffs

  Nblocks = floor(frames/(Nt*dec));
  printf("frames: %d Nblocks: %d\n", frames, Nblocks);

  unwrapped_dcts = zeros(Nblocks,Nt*K);
  rate_K_surface_ = zeros(frames, K);
  
  % create map on the fly from train database

  asurf = load("all_surf.txt"); [nr nc] = size(asurf);
  asurf = asurf(1:dec:nr,:);
  [map rms_map mx mx_ind unwrapped_dcts] = create_map_rms(asurf, Nt, K);
  %map = create_zigzag_map(Nt,K);

  %printf("non zero coeffs: %d\n", sum(sum(map == 1)));
  figure(2); clf;
  mesh(map);
  sumnz = zeros(1,Nblocks);
  dct2_sd = zeros(1,Nblocks);

  % create quantiser_num to r,c Luts

  rmap = cmap = zeros(1,Nt*K);
  for r=1:Nt
    for c=1:K
       quantiser_num = map(r,c);
       rmap(quantiser_num) = r;
       cmap(quantiser_num) = c;
    end
  end

  for n=1:Nblocks
    st = (n-1)*dec*Nt+1; en = st + dec*Nt - 1;
    %printf("st: %d en: %d\n", st, en);
    D = dct2(rate_K_surface(st:dec:en,:));

    % move over surface and work out quantiser
    % quantise, replace on map

    E = mapped = zeros(Nt,K);
    #{
    for r=1:Nt
      for c=1:K
        quantiser_num = map(r,c);
        if quantiser_num <= 40
          %printf("r %d c %d quantiser_num %d nlevels %d ", r, c, quantiser_num, nlevels(quantiser_num));
          %levels = quantisers(quantiser_num, 1:nlevels(quantiser_num));
          %quant_out = quantise(levels, D(r,c));
          E(r,c) = D(r,c);
          if E(r,c)
            mapped(r,c) = 1;
            sumnz(n)++;
          end
        end
      end
    end
    #}

    qn = 0;
    adct2_sd = mean(std(D-E));
    while adct2_sd > dist_dB
      qn++;
      E(rmap(qn), cmap(qn)) = 1*round(D(rmap(qn), cmap(qn))/1);
      adct2_sd = mean(std(D-E));
      %printf("qn %d %f\n", qn, adct2_sd);
    end
    sumnz(n) = qn;

    % note neat trick to interpolate to 10ms frames despite dec

    #{
    energy = sum(sum(E));
    Edc = E(1,1);
    E = E*1.2;
    E(1,1) = Edc/1.2;
    #}
    %E *= energy/sum(sum(E));
    rate_K_surface_(st:en,:) = idct2([sqrt(dec)*E; zeros(Nt*(dec-1), K)]);

    dct2_sd(n) = mean(std(D-E));
  end

  % figure(3); clf; mesh(mapped);
  figure(4); clf; plot(sumnz); hold on; plot([1 length(sumnz)],[mean(sumnz) mean(sumnz)]); hold off; title('Non Zero');
  figure(5); clf; plot(dct2_sd); title('DCT SD');
  printf("average dct spectral distortion: %3.2f dB\n", mean(dct2_sd));
  printf("mean number of coeffs/DCT: %3.2f/%d\n", mean(sumnz), Nt*K);
  printf("coeffs/second: %3.2f\n", mean(sumnz)/(Nt*Tf*dec));
  printf("bits/s: %3.2f\n", 2.9*mean(sumnz)/(Nt*Tf*dec));

  % optional 700C style post filter

  post_filter_en = 0;
  if post_filter_en
    for f=1:Nt*Nblocks
      mn = mean(rate_K_surface_(f,:));
      rate_K_surface_no_mean_(f,:) =  rate_K_surface_(f,:) - mn;
      rate_K_surface_(f,:) = mn + post_filter(rate_K_surface_no_mean_(f,:), rate_K_sample_freqs_kHz);
    end
  end

  % prevent /0 errors at end of run

  rate_K_surface_(dec*Nt*Nblocks+1:frames,:) = rate_K_surface(dec*Nt*Nblocks+1:frames,:); 
  model_ = resample_rate_L(model, rate_K_surface_, rate_K_sample_freqs_kHz);
  
  dist = std((rate_K_surface_(1:dec:frames,:) - rate_K_surface(1:dec:frames,:))');
  figure(1); clf; plot(dist); title('Rate K SD');
  printf("Rate K spectral distortion mean: %3.2f dB var: %3.2f\n", mean(dist), var(dist));
endfunction


% Basic unquantised rate K mel-sampling then back to rate L.  Now with "high end correction"

function [model_ rate_K_surface] = experiment_mel_freq(model, vq_en=0, plots=1, voicing)
  [frames nc] = size(model);
  K = 20; Fs = 16000; correct_rate_K_en = 0;

  for f=1:frames
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*Fs/(2000*pi);
    [rate_K_vec rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model(f,:), K, Fs);
    if correct_rate_K_en
      [tmp_ AmdB_] = resample_rate_L(model(f,:), rate_K_vec, rate_K_sample_freqs_kHz, Fs);
      [rate_K_vec_corrected orig_error error nasty_error_log nasty_error_m_log] = correct_rate_K_vec(rate_K_vec, rate_K_sample_freqs_kHz, AmdB, AmdB_, K, Wo, L, Fs);
      rate_K_surface(f,:) = rate_K_vec_corrected;
    else
      rate_K_surface(f,:) = rate_K_vec;
    end
  end

  model_ = resample_rate_L(model, rate_K_surface, rate_K_sample_freqs_kHz, Fs);
 
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


% Experiment with 10ms update rate for energy but 40ms update for spectrum,
% using linear interpolation for spectrum.

function model_c = experiment_energy_rate_linear(model, vq_en, plot_en)
  max_amp = 80;
  [frames nc] = size(model);

  % 10ms mel freq modelling and VQ

  model_ = experiment_mel_freq(model, vq_en, plot_en);

  % Remove energy.  Hmmmm, this is done on Ams rather than surface but that's
  % similar I guess

  e = zeros(1,frames);
  model_a = zeros(frames,max_amp+3);
  for f=1:frames
    L = min([model_(f,2) max_amp-1]);
    Am_ = model_(f,3:(L+2));
    AmdB_ = 20*log10(Am_);
    mean_f(f) = mean(AmdB_);
    AmdB_ -= mean_f(f);
    model_a(f,1) = model_(f,1); model_a(f,2) = L; model_a(f,3:(L+2)) = 10 .^ (AmdB_(1:L)/20);
  end

  % linear interp after removing energy (mean AmdB)

  model_b = experiment_dec_linear(model_a);

  % add back in energy

  model_c = zeros(frames,max_amp+3);
  for f=1:frames
    L = min([model_b(f,2) max_amp-1]);
    Am_ = model_b(f,3:(L+2));
    AmdB_ = 20*log10(Am_);
    AmdB_ += mean_f(f);
    model_c(f,1) = model_b(f,1); model_c(f,2) = L; model_c(f,3:(L+2)) = 10 .^ (AmdB_(1:L)/20);
  end

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



% Experimental AbyS decimator that chooses best frames to match
% surface based on AbyS approach.  Can apply post filter at different
% points, and optionally do fixed decimation, at rate K.  Didn't
% produce anything spectacular in AbyS mode, suggest another look with
% some sort of fbf display to see what's going on internally.
 
function model_ = experiment_dec_abys(model, M=8, vq_en=0, pf_en=1, fixed_dec=0, voicing)
  max_amp = 80;
  [frames nc] = size(model);
  model_ = zeros(frames, max_amp+3);

  printf("M: %d vq_en: %d pf_en: %d fixed_dec: %d\n", M, vq_en, pf_en, fixed_dec)

  % create frames x K surface

  K = 20; 
  [surface sample_freqs_kHz] = resample_const_rate_f_mel(model, K);
  target_surface = surface;

  % optionaly VQ surface

  if vq_en
    melvq;
    load train_120_vq; m=5;
       
    for f=1:frames
      mean_f(f) = mean(surface(f,:));
      rate_K_surface_no_mean(f,:) = surface(f,:) - mean_f(f);
    end

    [res rate_K_surface_ ind] = mbest(train_120_vq, rate_K_surface_no_mean, m);

    if pf_en == 1
      for f=1:frames
        rate_K_surface_(f,:) = post_filter(rate_K_surface_(f,:), sample_freqs_kHz, 1.5, voicing(f));
      end
    end
    
    for f=1:frames
      rate_K_surface_(f,:) += mean_f(f);
    end

    surface = rate_K_surface_;
  end

  % break into segments of M frames.  Fix end points, that two of the
  % frames we sample.  Then find best choice in between

  surface_ = zeros(frames, K);
  best_surface_ = zeros(frames, K);

  for f=1:M:frames-M

    left_vec = surface(f,:);
    right_vec = surface(f+M,:);
    resample_points = f:f+M-1;
    best_mse = 1E32;
    best_m = f+1;
    

    if fixed_dec
      m = f+M/2;
      printf("%d %d %d\n", f, m, M+f);
      centre_vec = surface(m,:);
      sample_points = [f m f+M];
      for k=1:K
        best_surface_(resample_points,k)  = interp1(sample_points, [left_vec(k) centre_vec(k) right_vec(k)], resample_points, "spline", 0);
      end 
    else
      printf("%d %d\n", f, M+f);
      for m=f+1:M+f-2

        % Use interpolation to construct candidate surface_ segment 
        % using just threee samples

        centre_vec = surface(m,:);
        sample_points = [f m f+M];
        mse = 0;
        for k=1:K
          surface_(resample_points,k)  = interp1(sample_points, [left_vec(k) centre_vec(k) right_vec(k)], resample_points, "spline", 0);
          mse += sum((target_surface(resample_points,k) - surface_(resample_points,k)).^2);
        end 

        % compare synthesised candidate to orginal and chose min ased on MSE

        if mse < best_mse
          best_mse = mse;
          best_m = m;
          best_surface_(resample_points,:) = surface_(resample_points,:);
        end

        printf("  m: %d mse: %f best_mse: %f best_m: %d\n", m, mse, best_mse, best_m);
      end
    end
  end

  if pf_en == 2

    % Optionally apply pf after interpolation, theory is interpolation
    % smooths spectrum in time so post filtering afterwards may be
    % useful.  Note we remove mean, this tends to move formats up and
    % anti-formants down when we multiply by a constant

    for f=1:frames
      mean_f = mean(best_surface_(f,:));
      rate_K_vec_no_mean = best_surface_(f,:) - mean_f;
      rate_K_vec_no_mean *= 1.2;
      best_surface_(f,:) = rate_K_vec_no_mean + mean_f;
    end
  end

  model_ = resample_rate_L(model, best_surface_, sample_freqs_kHz);
  figure(5);
  plot(mean_f,'+-')
  figure(6)
  hist(mean_f)
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


% smoothing using masking functions

function model_ = experiment_smoothed(model, bark_model=1)
  [frames nc] = size(model);

  model_ = model;

  for f=1:frames
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);
    [AmdB_ Am_freqs_kHz] = mask_model(AmdB, Wo, L, bark_model);
    model_(f,3:(L+2)) = 10 .^ ((AmdB_)/20);
  end

endfunction

#{
  
  My original idea was to used a 3-4 "resonators" to construct a
  piecewise model of the spectrum.  Kind of got distracted by the
  surface and mel sampling that ended up working OK.  This method was
  working OK, soem issues with background noise but rather easy to
  quantise.

  todo: get this working again
#}

function model_ = experiment_piecewise(model)

  [frames tmp] = size(model); max_amp = 160;

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


% helper function to plot surfaces

function plot_dft_surface(surface)
  [frames K] = size(surface);

  dft_surface = zeros(frames,K);
  for k=1:K
    dft_surface(:,k) = fft(surface(:,k).*hanning(frames));
  end

  % Set up a meaninful freq axis.  We only need real side. Sample rate
  % (frame rate) Fs=100Hz 

  Fs = 100; dF = Fs/frames;
  mesh(1:K, (1:frames/2)*dF, 20*log10(abs(dft_surface(1:frames/2,:,:))));
endfunction

