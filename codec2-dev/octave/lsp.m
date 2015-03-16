%
% Given input AR filter coefficients, this function calculates
% the corresponding Line spectral pair (LSP) parameters.
%
% Author : Madhukar Budagavi
%

1;

function w=atolsp(a)

  [m,n] = size(a);
  if m > n
    a = a';
  end
  m = length(a)-1;
  pz = [a,0] - [0 fliplr(a)];
  qz = [a,0] + [0 fliplr(a)];
  angpz = sort(angle(roots(pz)));
  angpz = angpz(find(angpz > 0));
  angqz = sort(angle(roots(qz)));
  angqz = angqz(find(angqz > 0 & angqz < pi));
  w = zeros(1,m);
  j = 1:(m/2);
  w(2*j) = angpz';
  w(2*j-1) = angqz';
endfunction

%
% Given input Line spectral pair (LSP) parameters of an AR filter,
% this function calculates the corresponding AR filter parameters
%
% Author : Madhukar Budagavi
%

function a=lsptoa(w)

  [m,n] = size(w);
  if m > n
    w = w';
  end
  m = length(w);
  pz = 1;
  qz = 1;
  for j = 1:(m/2),
    pz = conv(pz,[1 -2*cos(w(2*j)) 1]);
    qz = conv(qz,[1 -2*cos(w(2*j-1)) 1]);
  end
  pz = conv([1 -1],pz);
  qz = conv([1 1],qz);
  a = (pz+qz)/2;
  a = a(1:(m+1));

endfunction
