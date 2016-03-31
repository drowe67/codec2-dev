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

function dk_log = newamp_batch(samname, optional_Am_out_name, optional_Aw_out_name)
  newamp;
  more off;

  max_amp = 80;
  dec_in_freq = 1;
  postfilter = 1;
  dec_in_time = 0;
  synth_phase = 0;
  freq_quant = 0;
  amp_quant = 0;
  dk_log = [];

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);
  model_ = zeros(frames, nc);
  nom_masked_m = zeros(frames,max_amp);

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
   
  fam  = fopen(Am_out_name,"wb"); 
  if synth_phase
    faw = fopen(Aw_out_name,"wb"); 
  end

  % encoder loop ------------------------------------------------------

  %pp_bw = gen_pp_bw;

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
    a_non_masked_m = find(AmdB > maskdB);
    non_masked_m(f,1:length(a_non_masked_m)) = a_non_masked_m;

    if dec_in_freq
      [maskdB_ tmp1 D dk_] = decimate_in_freq(maskdB, 1);
      dk_log = [dk_log; dk_];
    end
    sd_sum += sum(maskdB - maskdB_);

    Am_ = zeros(1,max_amp);
    Am_ = 10 .^ (maskdB_(1:L-1)/20); 
    model_(f,3:(L+1)) = Am_;
  end
  printf("\nsd_sum: %5.2f\n", sd_sum/frames);

  % decoder loop -----------------------------------------------------

  for f=1:frames
    %printf("frame: %d\n", f);
    L = min([model_(f,2) max_amp-1]);
    Wo = model_(f,1);
    Am_ = model_(f,3:(L+2));

    maskdB_ = 20*log10(Am_);
    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;

    if dec_in_time
      % decimate mask samples in time
      maskdB_ = decimate_frame_rate(model_, 4, f, frames, mask_sample_freqs_kHz);

      % find turning points - prototype for finding PF freqs when we decimate in time

      d = maskdB_(2:L) - maskdB_(1:L-1);
      tp = [];
      for m=1:L-2
        if (d(m) > 0) && (d(m+1) < 0)
          tp = [tp m+1];
        end
      end
      a_non_masked_m = tp;
    else
      % read non-masked (PF freqs) from analysis stage
      % number of non-masked samples is variable when not using AbyS,
      % but fixed when using AbyS

      m = max(find(non_masked_m(f,:) > 0)); 
      a_non_masked_m = non_masked_m(f,1:m);
    end

    % post filter - bump up samples by 6dB, reduce mask by same level to normalise gain

    if postfilter

      % Apply post filter - enhances formants, suppresses others, as pe Part 1 blog
      % Pretty simple but makes a big difference

      maskdB_pf = maskdB_ - 6;
      maskdB_pf(a_non_masked_m) = maskdB_pf(a_non_masked_m) + 6;

      Am_ = zeros(1,max_amp);
      Am_ = 10 .^ (maskdB_pf(1:L-1)/20); 
      model_(f,3:(L+1)) = Am_;
    else
      maskdB_pf = maskdB_;
    end

    Am_ = zeros(1,max_amp);
    Am_(2:L) = 10 .^ (maskdB_pf(1:L-1)/20);  % C array doesnt use A[0]
    fwrite(fam, Am_, "float32");

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
  if synth_phase
    fclose(faw);
  end

  printf("\n");

endfunction


