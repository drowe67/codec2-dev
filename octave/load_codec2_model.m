% load_codec2_model.m
%
% David Rowe 2023
%
% Loads binary model files from c2sim:
%  $ ./c2sim ../raw/big_dog.raw --modelout big_dog_model.bin
%  octave> model=load_codec2_model('big_dog_model.bin');
%
% model retuned is compatable with "c2sim --dump" sample_model.txt dump file

function model = load_codec2_model(fn)
  max_amp = 160; f = 1;
  
  s = dir(fn);         
  filesize = s.bytes;
  nrows = s.bytes/1300;
  ncols = (1 + 1 + (max_amp/2) + (max_amp/2) + 1);
  model = zeros(nrows, ncols);
       
  fp = fopen(fn,"rb");
  do
    [Wo count] = fread(fp, 1, "float");
    if count
      L = fread(fp, 1, "int");
      A = fread(fp, max_amp+1, "float")';
      phi = fread(fp, max_amp+1, "float")';
      voiced = fread(fp, 1, "int");
      amodel = [Wo L A(2:end) voiced];
      model(f,:) = amodel;
      f++;
    end  
  until (count == 0)
  close f;
endfunction
