% newamp_batch.m
%
% Copyright David Rowe 2015
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Octave script to batch process model parameters using the new
% amplitude model.  Used for generating samples we can listen to.
%
% Usage:
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
%   $ cd ~/codec2-dev/octave
%   octave:14> newamp_batch("../build_linux/src/hts1a")
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --amread hts1a_am.out -o - | play -t raw -r 8000 -s -2 -


% process a whole file and write results

function [dk_log D1_log] = newamp_batch(samname, optional_Am_out_name, optional_Aw_out_name)
  newamp;
  more off;

  max_amp = 80;
  dec_in_freq = 1;
  postfilter = 0;
  dec_in_time = 1;
  synth_phase = 1;
  vq_en = 1;
  D_log = []; dk_log = []; D1_log = []; ind_log = [];
  train = 0;

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);
  model_ = zeros(frames, nc);
  non_masked_m = zeros(frames,max_amp);

  if nargin == 1
    Am_out_name = sprintf("%s_am.out", samname);
    Aw_out_name = sprintf("%s_aw.out", samname);
  end
  if nargin >= 2
    Am_out_name = optional_Am_out_name;
  end
  if nargin >= 3
    Aw_out_name = optional_Aw_out_name;
    faw = fopen(Aw_out_name,"wb"); 
  end
   
  if vq_en
    load vq;
  end


  % encoder loop ------------------------------------------------------

  sd_sum = 0;
  for f=1:frames
    printf("%d ", f);   
    model_(f,2) = L = min([model(f,2) max_amp-1]);
    model_(f,1) = Wo = model(f,1);
    model_(f,3:(L+2)) = Am = model(f,3:(L+2));

    AmdB = 20*log10(Am);

    % find mask

    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    maskdB = mask_model(AmdB, Wo, L);

    % save post filter positions for decode

    maskdB_ = maskdB;
    if (postfilter == 1) || (postfilter == 3)
      a_non_masked_m = find(AmdB > maskdB);
    end
    if postfilter == 2
      a_non_masked_m = est_pf_locations(maskdB_);
    end
    if postfilter
      if length(a_non_masked_m)
        non_masked_m(f,1:length(a_non_masked_m)) = a_non_masked_m;
      end
    end

    if postfilter == 3
       maskdB_ = maskdB_ - 9;
       maskdB_(a_non_masked_m) = maskdB_(a_non_masked_m) + 9;
    end

    AmdB_ = AmdB;
    if dec_in_freq
      if vq_en
        [AmdB_ tmp1 Dabs dk_ D1 ind] = decimate_in_freq(maskdB, 1, 10, vq);
        ind_log = [ind_log; ind];
      else
        [AmdB_ tmp1 D dk_ D1] = decimate_in_freq(maskdB, 0, 10);
      end
      if train
        dk_log = [dk_log; dk_];
        D1_log = [D1_log D1];
      end
    end
    sd_sum += std(maskdB - AmdB_);

    Am_ = zeros(1,max_amp);
    Am_ = 10 .^ (AmdB_(1:L-1)/20); 
    model_(f,3:(L+1)) = Am_;
  end

  if train == 0
    decode_model(model_, Am_out_name, Aw_out_name, synth_phase, dec_in_time);
  end

  printf("\nsd_sum: %5.2f\n", sd_sum/frames);
  printf("\n");
endfunction


% generate array of indexes
% convert to bits and save to file of one char/bit
% function to encode and save to file of bits
% function to decode from file of bits
% quantise Wo
% move voicing bit through, move Wo through, maybe just load Wo and L and v?

function decode_model(model_, Am_out_name, Aw_out_name, synth_phase, dec_in_time)
  max_amp = 80;

  fam  = fopen(Am_out_name,"wb"); 
  if synth_phase
    faw = fopen(Aw_out_name,"wb"); 
  end

  % decoder loop -----------------------------------------------------

  [frames tmp] = size(model_);
  for f=1:frames
    %printf("frame: %d\n", f);
    L = min([model_(f,2) max_amp-1]);
    Wo = model_(f,1);
    Am_ = model_(f,3:(L+2));
    AmdB_ = 20*log10(Am_);
    sample_freqs_kHz = (1:L)*Wo*4/pi;

    % run post filter ahead of time so dec in time has post filtered frames to work with

    if f+3 <= frames
      model_(f+3,:) = post_filter(model_(f+3,:));
    end

    if dec_in_time
      % decimate mask samples in time

      decimate = 4;
      AmdB_ = decimate_frame_rate(model_, decimate, f, frames, sample_freqs_kHz);
    end

    Am_ = zeros(1,max_amp);
    Am_(2:L) = 10 .^ (AmdB_(1:L-1)/20);  % C array doesnt use A[0]
    fwrite(fam, Am_, "float32");

    if synth_phase

      % synthesis phase spectra from magnitiude spectra using minimum phase techniques

      fft_enc = 512;
      model_(f,3:(L+2)) = 10 .^ (AmdB_(1:L)/20);
      phase = determine_phase(model_, f);
      assert(length(phase) == fft_enc);
      Aw = zeros(1, fft_enc*2); 
      Aw(1:2:fft_enc*2) = cos(phase);
      Aw(2:2:fft_enc*2) = -sin(phase);

      fwrite(faw, Aw, "float32");    
    end
  end

  fclose(fam);
  fclose(faw);
endfunction


