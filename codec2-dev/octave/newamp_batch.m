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
  fully_quant = 1;
  Wo_quant = 1;

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

    if Wo_quant || fully_quant
      ind_Wo = encode_log_Wo(Wo, 6);    
      model_(f,1) = Wo_ = decode_log_Wo(ind_Wo, 6);
      L_  = floor(pi/Wo);
      model_(f,2) = L_ = min([L_ max_amp-1]);
    else
      model_(f,1) = Wo_ = Wo;
      model_(f,2) = L_ = L;
    end

    % find mask

    mask_sample_freqs_kHz = (1:L_)*Wo_*4/pi;
    maskdB = mask_model(AmdB, Wo_, L_);
    maskdB_ = maskdB;
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

    if fully_quant
      if vq_en
        ind_log = [ind_log; ind_Wo voicing(f) (ind_vq-1) 0];
      end
    end

    Am_ = zeros(1,max_amp);
    Am_ = 10 .^ (AmdB_(1:L_)/20); 
    model_(f,3:(L_+2)) = Am_;
  end
  
  if train == 0
    if fully_quant
      bit_stream_name = strcat(samname,".bit");
      bits_per_param = [6 1 8 8 4 1];
      write_bit_stream_file(bit_stream_name, ind_log, bits_per_param);
      decode_bit_stream_file(samname, bits_per_param);
    else
      decode_model(model_, samname, synth_phase, dec_in_time);
    end
  end

  printf("\nsd_sum: %5.2f\n", sd_sum/frames);
  printf("\n");
endfunction

