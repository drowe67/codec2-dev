% tnewamp1.m
%
% Copyright David Rowe 2017
% This program is distributed under the terms of the GNU General Public License 
% Version 2

#{

  Octave script to compare Octave and C versions of newamp1 processing, in order to test C port.

  c2sim -> dump files -> $ ../build_linux/unittest/tnewamp1 -> octave:1> tnewamp1
  Usage:

    1/ build codec2 with -DDUMP - see codec2-dev/README

    2/ Generate dump files using c2sim (just need to do this once)
       $ cd codec2-dev/build_linux/src
       $ ./c2sim ../../raw/hts1a.raw --phase0 --postfilter --dump hts1a --lpc 10 --dump_pitch_e hts1a_pitche.txt

    3/ Run C version which generates a file of Octave test vectors as ouput:

      $ cd codec2-dev/build_linux/unittest
      $ ./tnewamp1 ../../raw/hts1a
            
    4/ Run Octave script to generate Octave test vectors and compare with C.

      octave:1> tnewamp1("../build_linux/src/hts1a")

#}


% In general, this function processes a bunch of amplitudes, we then
% use c2sim to hear the results.  Bunch of different experiments below

function tnewamp1(input_prefix)
  newamp;
  autotest;
  more off;

  max_amp = 80;
  postfilter = 0;   % optional postfiler that runs on Am, not used atm
  synth_phase = 1;

  if nargin == 1
    output_prefix = input_prefix;
  end
  model_name = strcat(input_prefix,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);

  voicing_name = strcat(input_prefix,"_pitche.txt");
  voicing = zeros(1,frames);
  
  if exist(voicing_name, "file") == 2
    pitche = load(voicing_name);
    voicing = pitche(:, 3);
  end

  % Load in C vectors and compare -----------------------------------------
 
  load("../build_linux/unittest/tnewamp1_out.txt");
  
  K = 20;
  [frames tmp] = size(rate_K_surface_c);
  [rate_K_surface sample_freqs_kHz] = resample_const_rate_f_mel(model(1:frames,:), K);

  melvq;
  load train_120_vq; m=5;
       
  for f=1:frames
    mean_f(f) = mean(rate_K_surface(f,:));
    rate_K_surface_no_mean(f,:) = rate_K_surface(f,:) - mean_f(f);
  end

  [res rate_K_surface_no_mean_ ind] = mbest(train_120_vq, rate_K_surface_no_mean, m);

  for f=1:frames
    rate_K_surface_no_mean_(f,:) = post_filter(rate_K_surface_no_mean_(f,:), sample_freqs_kHz, 1.5);
  end
    
  rate_K_surface_ = zeros(frames, K);
  interpolated_surface_ = zeros(frames, K);
  energy_q = create_energy_q;
  M = 4;
  for f=1:frames   
    [mean_f_ indx] = quantise(energy_q, mean_f(f));
    indexes(f,3) = indx - 1;
    rate_K_surface_(f,:) = rate_K_surface_no_mean_(f,:) + mean_f_;
  end

  % simulated decoder
  % break into segments of M frames.  We have 2 samples spaced M apart
  % and interpolate the rest.

  Nfft_phase = 128;
  model_ = zeros(frames, max_amp+2);
  voicing_ = zeros(1,frames);
  H = zeros(frames, max_amp);
  model_(1,1) = Wo_left = 2*pi/100;
  voicing_left = 0;
  left_vec = zeros(1,K);

  for f=1:M:frames   

    if voicing(f)
      index = encode_log_Wo(model(f,1), 6);
      if index == 0
        index = 1;
      end
      model_(f,1) = decode_log_Wo(index, 6);
    else
      model_(f,1) = 2*pi/100;
    end

    Wo_right = model_(f,1);
    voicing_right = voicing(f);
    [Wo_ avoicing_] = interp_Wo_v(Wo_left, Wo_right, voicing_left, voicing_right);

    if f > M
      model_(f-M:f-1,1) = Wo_;
      voicing_(f-M:f-1) = avoicing_;
      model_(f-M:f-1,2) = floor(pi ./ model_(f-M:f-1,1)); % calculate L for each interpolated Wo
    end

    right_vec = rate_K_surface_(f,:);

    if f > M
      sample_points = [f-M f];
      resample_points = f-M:f-1;
      for k=1:K
        interpolated_surface_(resample_points,k) = interp_linear(sample_points, [left_vec(k) right_vec(k)], resample_points);
      end

      for k=f-M:f-1
        model_(k,:) = resample_rate_L(model_(k,:), interpolated_surface_(k,:), sample_freqs_kHz);
        phase = determine_phase(model_, k, Nfft_phase);
        for m=1:model_(k,2)
          b = round(m*model_(k,1)*Nfft_phase/(2*pi));  % map harmonic centre to DFT bin
          H(k,m) = exp(-j*phase(b+1));
        end     
      end
   end
   
   % update for next time

   Wo_left = Wo_right;
   voicing_left = voicing_right;
   left_vec = right_vec;
   
  end

  %model_(1,1:77)
  %model__c(1,1:77)
  %sum(model_(1,1:77)-model__c(1,1:77))
  %[mx mxi] = max(model_(1,1:77)-model__c(1,1:77))

  %interpolated_surface_(1,:)
  %interpolated_surface__c(1,:)
  %sum(interpolated_surface_(1,:) - interpolated_surface__c(1,:))


  %Hm(2,:) - Hm_c(2,:)
  
  figure(1);
  mesh(angle(H));
  figure(2);
  mesh(angle(H_c));
  figure(3);
  mesh(abs(H - H_c));

  check(rate_K_surface, rate_K_surface_c, 'rate_K_surface', 0.01);
  check(mean_f, mean_c, 'mean', 0.01);
  check(rate_K_surface_, rate_K_surface__c, 'rate_K_surface_', 0.01);
  check(interpolated_surface_, interpolated_surface__c, 'interpolated_surface_', 0.01);
  check(model_(:,1), model__c(:,1), 'interpolated Wo_', 0.001);
  check(voicing_, voicing__c, 'interpolated voicing');
  check(model_(:,3:max_amp+2), model__c(:,3:max_amp+2), 'rate L Am surface ', 0.1);
  check(H, H_c, 'phase surface');

  #{

  for f=1:frames
    printf("%d ", f);   
    Wo = model_(f,1);
    L = min([model_(f,2) max_amp-1]);
    Am = model_(f,3:(L+2));

    Am_ = zeros(1,max_amp);
    Am_(2:L) = Am(1:L-1);

    % optional post filter on linear {Am}, boosts higher amplitudes more than lower,
    % improving shape of formants and reducing muffling.  Note energy
    % normalisation

    if postfilter
      e1 = sum(Am_(2:L).^2);
      Am_(2:L) = Am_(2:L) .^ 1.5;         
      e2 = sum(Am_(2:L).^2);
      Am_(2:L) *= sqrt(e1/e2);
    end

    fwrite(fam, Am_, "float32");
    fwrite(fWo, Wo, "float32");

    if synth_phase

      % synthesis phase spectra from magnitiude spectra using minimum phase techniques

      fft_enc = 128;
      phase = determine_phase(model_, f, fft_enc);
      assert(length(phase) == fft_enc);
      %Aw = zeros(1, fft_enc*2); 
      %Aw(1:2:fft_enc*2) = cos(phase);
      %Aw(2:2:fft_enc*2) = -sin(phase);

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

  #}

  printf("\n")

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
  energy_q = 10 + 40/16*(0:15);
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
      voicing(f) = 1;
      model_(f,1) = decode_log_Wo(indexes(f,4), 6);
    else
      indexes(f,4) = 0;
      voicing(f) = 0;
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
