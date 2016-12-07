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

function [fvec_log amps_log] = newamp1_batch(samname, optional_Am_out_name, optional_Aw_out_name)
  newamp;
  more off;

  max_amp = 80;
  load vq;

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);

  if nargin == 2
    Am_out_name = optional_Am_out_name;
  else
    Am_out_name = sprintf("%s_am.out", samname);
  end

  fam  = fopen(Am_out_name,"wb"); 

  % encoder loop ------------------------------------------------------

  fvec_log = []; amps_log = [];

  for f=1:frames
    printf("%d ", f);   

    Wo = model(f,1);
    L = min([model(f,2) max_amp-1]);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);

    % find mask

    #{
    maskdB = mask_model(AmdB, Wo, L);
    AmdB_ = maskdB;
    [mx mx_ind] = max(AmdB_);
    AmdB_(mx_ind) += 6;
    #}

    [AmdB_ res fvec fvec_ amps] = piecewise_model(AmdB, Wo, vq, 1);
    fvec_log = [fvec_log; fvec];
    amps_log = [amps_log; amps];
    #{
    l1000 = floor(L/4);     
    AmdB_ = AmdB;
    [mx mx_ind] = max(AmdB_(1:l1000));
    mask_sample_freqs_kHz = (1:l1000)*Wo*4/pi;
    AmdB_(1:l1000) = parabolic_resonator(mx_ind*Wo*4/pi, mask_sample_freqs_kHz) + mx;  
    #}

    Am_ = zeros(1,max_amp);
    Am_(2:L) = 10 .^ (AmdB_(1:L-1)/20);  % C array doesnt use A[0]
    fwrite(fam, Am_, "float32");
    Am_ = zeros(1,max_amp);
  end

  fclose(fam);
  printf("\n")

endfunction
  
