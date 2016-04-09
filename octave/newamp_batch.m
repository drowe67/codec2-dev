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
% 

% Generate voicing file:
% $ ./c2sim ../../raw/hts2a.raw --dump hts2a --lpc 10 --phase0 --dump_pitch_e hts2a_pitche.txt

% Bit stream decode:
% $ ./c2sim ../../raw/hts2a.raw --amead hts2a_am.out --awread hts2a_aw.out --Woread hts2a_Wo.out --phase0 --postfilter -o - | play -t raw -r 8000 -s -2 -


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
  %frames = 5;

  voicing_name = strcat(samname,"_pitche.txt");
  if exist(voicing_name, "file") == 2
    pitche = load(voicing_name);
    voicing = pitche(:, 3);
  end
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
    Wo = model(f,1); L = model(f,2);
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
        [AmdB_ tmp1 Dabs dk_ D1 ind_vq] = decimate_in_freq(maskdB, 1, 10, vq);
      else
        [AmdB_ tmp1 D dk_ D1] = decimate_in_freq(maskdB, 0, 10);
      end
      if train
        dk_log = [dk_log; dk_];
        D1_log = [D1_log D1];
      end
    end
    sd_sum += std(maskdB - AmdB_);

    ind_Wo = encode_log_Wo(Wo, 6);    
    if vq_en
      ind_log = [ind_log; ind_Wo voicing(f) (ind_vq-1) 0];
    end
    model_(f,1) = Wo_ = decode_log_Wo(ind_Wo, 6);
    L_  = floor(pi/Wo);
    model_(f,2) = L_ = min([L_ max_amp-1]);
    L_ = min(L_, length(AmdB_));
    Am_ = zeros(1,max_amp);
    Am_ = 10 .^ (AmdB_(1:L_)/20); 
    model_(f,3:(L_+2)) = Am_;
    printf("Wo: %f Wo_: %f L: %d L_: %d\n", Wo, Wo_, L, L_);
  end

  if train == 0
    bit_stream_name = strcat(samname,".bit");
    bits_per_param = [6 1 8 8 4 1];
    write_bit_stream_file(bit_stream_name, ind_log, bits_per_param);
    %decode_model(model_, samname, synth_phase, dec_in_time);
    decode_bit_stream_file(samname, bits_per_param);
  end

  printf("\nsd_sum: %5.2f\n", sd_sum/frames);
  printf("\n");
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
      printf("i: %d bits_per_param: %d\n", i, bits_per_param(i));
      some_bits = index_to_bits(arow(i), bits_per_param(i));
      frame_of_bits = [frame_of_bits some_bits];
    end
    fwrite(fbit, frame_of_bits, "uchar");
    arow
    frame_of_bits

  end

  fclose(fbit);
endfunction


% Decode from a bit stream file

function decode_bit_stream_file(samname, bits_per_param)
  max_amp = 80;
  nc = max_amp + 3;
  load vq;
  [tmp1 k2 tmp2] = size(vq)
  k = k2/2;

  bit_stream_name = strcat(samname,".bit")
  fbit  = fopen(bit_stream_name, "rb"); 

  model_= []; log_v = [];
  nind = length(bits_per_param);
  bits_per_frame = sum(bits_per_param);

  % read a frame, decode to indexes, fill in model_ array

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
    
    % convert index to model parameters

    amodel_ = zeros(1,nc);
    amodel_(1) = Wo_ = decode_log_Wo(ind(1), 6);
    L_  = floor(pi/Wo_);
    amodel_(2) = L_ = min([L_ max_amp-1]);

    [dk_ D1_] = index_to_params(ind(3:5)+1, vq);
    AmdB_ = params_to_mask(L_, k, dk_, D1_);

    Am_ = zeros(1,max_amp);
    Am_ = 10 .^ (AmdB_(1:L_)/20); 
    amodel_(3:(L_+2)) = Am_;
    model_ = [ model_; amodel_; zeros(3,nc)];

    log_v = [log_v ind(2)];

    % read next frame

    [frame nread] = fread(fbit, sum(bits_per_param), "uchar");
  endwhile

  % decode entire array of model parameters

  decode_model(model_, samname, 1, 1);

  % save voicing file
  
  v_out_name = sprintf("%s_v.txt", samname);
  fv  = fopen(v_out_name,"wt"); 
  for f=1:length(log_v)
    for i=1:4
      fprintf(fv,"%d\n",log_v(f));
    end
  end
  fclose(fv);
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


function decode_model(model_, samname, synth_phase, dec_in_time)
  max_amp = 80;

  Am_out_name = sprintf("%s_am.out", samname);
  Aw_out_name = sprintf("%s_aw.out", samname);
  Wo_out_name = sprintf("%s_Wo.out", samname);
  fam  = fopen(Am_out_name,"wb"); 
  fWo  = fopen(Wo_out_name,"wb"); 
  if synth_phase
    faw = fopen(Aw_out_name,"wb"); 
  end

  % decoder loop -----------------------------------------------------

  [frames tmp] = size(model_)
  
  for f=1:frames
    %printf("frame: %d\n", f);
    L = min([model_(f,2) max_amp-1]);
    Wo = model_(f,1);
    Am_ = model_(f,3:(L+2));
    AmdB_ = 20*log10(Am_);
    sample_freqs_kHz = (1:L)*Wo*4/pi;
    %printf("Wo: %f Wo_: %f L: %d L_: %d\n", Wo, Wo_, L, L_);

    % run post filter ahead of time so dec in time has post filtered frames to work with

    if f+3 <= frames
      model_(f+3,:) = post_filter(model_(f+3,:));
    end

    if dec_in_time
      % decimate mask samples in time

      decimate = 4;
      [AmdB_ Wo L] = decimate_frame_rate(model_, decimate, f, frames, sample_freqs_kHz);
    end

    Am_ = zeros(1,max_amp);
    Am_(2:L) = 10 .^ (AmdB_(1:L-1)/20);  % C array doesnt use A[0]
    fwrite(fam, Am_, "float32");
    fwrite(fWo, Wo, "float32");

    if synth_phase

      % synthesis phase spectra from magnitiude spectra using minimum phase techniques

      fft_enc = 512;
      model_(f,1) = Wo;
      model_(f,2) = L;
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
  fclose(fWo);
  fclose(faw);
endfunction


