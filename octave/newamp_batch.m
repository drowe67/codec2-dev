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

function newamp_batch(samname, optional_Am_out_name, optional_Aw_out_name)
  newamp;
  more off;

  max_amp = 80;
  decimate_in_freq = 1;
  postfilter = 1;
  decimate_in_time = 0;
  synth_phase = 0;

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);

  if nargin == 1
    Am_out_name = sprintf("%s_am.out", samname);
    Aw_out_name = sprintf("%s_aw.out", samname);
  end
  if nargin >= 2
    Am_out_name = optional_Am_out_name;
  end
  if nargin >= 3
    Aw_out_name = optional_Aw_out_name;
  end
   
  fam  = fopen(Am_out_name,"wb"); 
  faw = fopen(Aw_out_name,"wb"); 

  for f=1:frames
    printf("%d ", f);
    L = min([model(f,2) max_amp-1]);
    Wo = model(f,1);
    Am = model(f,3:(L+2));

    AmdB = 20*log10(Am);

    % find mask

    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    maskdB = mask_model(AmdB, Wo, L);

    if decimate_in_freq
      % decimate in mask samples in frequency
      [decmaskdB decmasksamples] = make_decmask_abys(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz);
      non_masked_m = decmasksamples(:,2);
      maskdB_ = decmaskdB;
    else
      % just approximate decimation in freq by using those mask samples we can hear
      maskdB_ = maskdB;
      non_masked_m = find(AmdB > maskdB);
    end

    if decimate_in_time
      % decimate mask samples in time
      maskdB_ = decimate_frame_rate(maskdB_, model, decimate, f, frames, mask_sample_freqs_kHz);
    end

    % post filter - bump up samples by 6dB, reduce mask by same level to normalise gain

    if postfilter
      maskdB_pf = maskdB_ - 6;
      maskdB_pf(non_masked_m) = maskdB_pf(non_masked_m) + 6;
    else
      maskdB_pf = maskdB_;
    end

    if 0 
      % Early work as per blog post part 1
      % Attempt 1
      maskdB_pf = zeros(1,L);
      maskdB_pf(non_masked_m) = maskdB(non_masked_m);
      % Attempt 2
      %maskdB_pf = maskdB;
      % Attempt 3
      %maskdB_pf = maskdB;
      %maskdB_pf(non_masked_m) += 6;
    end

    Am_ = zeros(1,max_amp);
    Am_(2:L) = 10 .^ (maskdB_pf(1:L-1)/20);  % C array doesnt use A[0]
    fwrite(fam, Am_, "float32");

    if synth_phase
      % synthesis phase spectra from magnitiude spectra using minimum phase techniques
      fft_enc = 512;
      phase = determine_phase(model, f);
      assert(length(phase) == fft_enc);
      Aw = zeros(1, fft_enc*2); 
      Aw(1:2:fft_enc*2) = cos(phase);
      Aw(2:2:fft_enc*2) = -sin(phase);
 
      fwrite(faw, Aw, "float32");    
    end
  end

  fclose(fam);
  fclose(faw);

  printf("\n");

endfunction


