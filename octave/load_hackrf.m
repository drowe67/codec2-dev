% load_hackrf.m
%
% David Rowe Oct 2015

function s = load_hackrf(fn, n)
  fs = fopen(fn,"rb");
  if nargin == 1
    iq = fread(fs,Inf,"schar");
  else
    iq = fread(fs,2*n,"schar");
  end  
  fclose(fs);
  l = length(iq);
  s =  iq(1:2:l) + j*iq(2:2:l);
endfunction
