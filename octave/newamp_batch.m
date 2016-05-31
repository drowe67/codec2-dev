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
  dec_in_time = 1;
  dec_in_freq = 1;
  decimate = 4;
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
      dk_log = [dk_log; dk_];
      D1_log = [D1_log; D1_];
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
  decode_from_bit_stream(samname);

endfunction
  
