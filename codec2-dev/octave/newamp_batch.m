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

function newamp_batch(samname)
  newamp;
  
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);
  max_amp = 80;

  Am_out_name = sprintf("%s_am.out", samname);
  fam = fopen(Am_out_name,"wb"); 

  for f=1:frames
    
    L = min([model(f,2) max_amp-1]);
    Wo = model(f,1);
    Am = model(f,3:(L+2));
    AmdB = 20*log10(Am);

    maskdB = mask_model(AmdB, Wo, L);
    mask_sample_freqs_kHz = (1:L)*Wo*4/pi;
    [newmaskdB local_maxima] = make_newmask(maskdB, AmdB, Wo, L, mask_sample_freqs_kHz);

    [nlm tmp] = size(local_maxima);
    non_masked_m = local_maxima(1:min(4,nlm),2);
    maskdB_pf = newmaskdB;
    maskdB_pf(non_masked_m) = maskdB_pf(non_masked_m) + 6;
    
    Am_ = zeros(1,max_amp);
    Am_(2:L) = 10 .^ (maskdB_pf(1:L-1)/20);
    
    % save to file
    fwrite(fam, Am_, "float32");
  end

  fclose(fam);

endfunction

