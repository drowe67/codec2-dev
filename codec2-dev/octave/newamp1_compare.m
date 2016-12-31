% newamp1_compare.m
%
% Copyright David Rowe 2016
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Compare model, Wo, and voicing files, used for checking refactoring and C port


function newamp1_compare(prefixa, prefixb)
  autotest;

  frames = 100;
  a = load_params(prefixa, frames);
  b = load_params(prefixb, frames);
 
  check(a.Am, b.Am, 'Am');
  check(a.Hm, b.Hm, 'Hm');
  check(a.Wo, b.Wo, 'Wo');
  check(a.v, b.v, 'v');

  figure(1); clf; plot(a.v); hold on; plot(a.v-b.v,'r'); hold off;
endfunction


function params = load_params(prefix, frames)
  max_amp = 80;
  fft_enc = 512;

  Am_out_name = sprintf("%s_am.out", prefix);
  fam  = fopen(Am_out_name,"rb");
  Wo_out_name = sprintf("%s_Wo.out", prefix);
  fWo  = fopen(Wo_out_name,"rb"); 
  Hm_out_name = sprintf("%s_hm.out", prefix);
  fhm = fopen(Hm_out_name,"rb"); 

  % load up values from binary files

  params.Am = zeros(frames, max_amp);
  params.Wo = zeros(frames, 1);
  params.v = zeros(frames, 1);
  params.Hm = zeros(frames, max_amp);
  for f=1:frames
    params.Am(f,:) = fread(fam, max_amp, "float32");
    params.Wo(f)   = fread(fWo, 1, "float32");
    params.Hm(f,:) = fread(fhm, max_amp, "float32");
  end

  fclose(fam); fclose(fWo); fclose(fhm);

  % voicing is a text file

  v_out_name = sprintf("%s_v.txt", prefix);
  params.v = load(v_out_name);

endfunction
