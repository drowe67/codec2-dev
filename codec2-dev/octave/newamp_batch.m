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

function newamp_batch(samname, optional_Am_out_name)
  newamp;
  more off;

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);
  max_amp = 80;

  if nargin == 1
    Am_out_name = sprintf("%s_am.out", samname);
  end
  if nargin == 2
    Am_out_name = optional_Am_out_name;
  end
    
  fam = fopen(Am_out_name,"wb"); 

  for f=1:frames
    
    L = min([model(f,2) max_amp-1]);
    Wo = model(f,1);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);

    maskdB = mask_model(AmdB, Wo, L);
    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    [newmaskdB local_maxima] = make_decmask_abys(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz);

    [nlm tmp] = size(local_maxima);
    non_masked_m = local_maxima(1:min(4,nlm),2);

    % post filter - bump up samples by 6dB, reduce mask by same level to normalise gain

    maskdB_pf = newmaskdB - 6;
    maskdB_pf(non_masked_m) = maskdB_pf(non_masked_m) + 6;

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
    Am_(2:L) = 10 .^ (maskdB_pf(1:L-1)/20);
    
    % save to file
    fwrite(fam, Am_, "float32");
  end

  fclose(fam);

endfunction

