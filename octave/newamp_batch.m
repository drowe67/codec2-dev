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
% Or with a little more processing:
%   codec2-dev/build_linux/src$ ./c2sim ../../raw/hts2a.raw --amread hts2a_am.out --awread hts2a_aw.out --phase0 --postfilter --Woread hts2a_Wo.out -o - | play -q -t raw -r 8000 -s -2 -


% process a whole file and write results

function [dk_log D1_log] = newamp_batch(samname, optional_Am_out_name, optional_Aw_out_name)
  newamp;
  more off;

  max_amp = 80;
  k = 10;
  decimate = 4;
  dec_in_freq = 1;
  dec_in_time = 1;
  synth_phase = 1;
  vq_en = 1;
  dk_log = []; D1_log = [];
  train = 0;
  Wo_quant = 1;
  ind_log = [];

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
   
  fam  = fopen(Am_out_name,"wb"); 
  if synth_phase
    faw = fopen(Aw_out_name,"wb"); 
  end

  Wo_out_name = sprintf("%s_Wo.out", samname);
  fWo  = fopen(Wo_out_name,"wb"); 

  voicing_name = strcat(samname,"_pitche.txt");
  voicing = zeros(1,frames);
  if exist(voicing_name, "file") == 2
    pitche = load(voicing_name);
    voicing = pitche(:, 3);
  end

  if vq_en
    load vq;
  end

  % encoder loop ------------------------------------------------------

  sd_sum = 0;
  for f=1:frames
    printf("%d ", f);   

    Wo = model(f,1);
    L = min([model(f,2) max_amp-1]);

    if Wo_quant
      ind_Wo = encode_log_Wo(Wo, 6);    
      Wo = decode_log_Wo(ind_Wo, 6);
      L = floor(pi/Wo);
      L = min([L max_amp-1]);
    end

    model_(f,1) = Wo;
    model_(f,2) = L;
    model_(f,3:(L+2)) = Am = model(f,3:(L+2));

    AmdB = 20*log10(Am);

    % find mask

    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    maskdB = mask_model(AmdB, Wo, L);

    maskdB_ = maskdB;
    if dec_in_freq
      if vq_en
        [maskdB_ tmp1 D dk_ D1_ ind_vq] = decimate_in_freq(maskdB, 1, k, vq);
      else
        [maskdB_ tmp1 D dk_ D1_] = decimate_in_freq(maskdB, 1, k);
      end
      if train
        dk_log = [dk_log; dk_];
        D1_log = [D1_log; D1_];
      end
    end
    %maskdB_pf = maskdB_*1.5;
    %maskdB_pf += max(maskdB_) - max(maskdB_pf);

    % log info for bit stream

    ind_log = [ind_log; ind_Wo voicing(f) (ind_vq-1) 0];

    %sd_sum += sum(maskdB_ - maskdB);

    Am_ = zeros(1,max_amp);
    Am_ = 10 .^ (maskdB_(1:L)/20); 
    model_(f,3:(L+2)) = Am_;
  end

  bit_stream_name = strcat(samname,".bit");
  bits_per_param = [6 1 8 8 4 1];
  write_bit_stream_file(bit_stream_name, ind_log, bits_per_param);

  % read in bit stream and convert to ind_log[]

  ind_log = [];
  fbit  = fopen(bit_stream_name, "rb"); 
  bits_per_frame = sum(bits_per_param);
  nind = length(bits_per_param);

  [frame nread] = fread(fbit, sum(bits_per_param), "uchar");
  while (nread == bits_per_frame)
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

  % convert ind_log to modem params

  frames = 4*length(ind_log);
  model_ = zeros(frames, max_amp+2);
  v      = zeros(frames,1);

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

  if exist(voicing_name, "file") == 2
    % save voicing file
  
    v_out_name = sprintf("%s_v.txt", samname);
    fv  = fopen(v_out_name,"wt"); 
    for f=1:length(v)
      fprintf(fv,"%d\n", v(f));
    end
    fclose(fv);
  end

  printf("\nsd_sum: %5.2f\n", sd_sum/frames);
  printf("\n");
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
