% gen_complex_short.m
%
% David Rowe Feb 2017
%
% Generates a complex short signal for HackRF testing.

Fs = 8E3;
T = 10;
A = 32000;
f = 1000;

N = T*Fs;
t = 0:N-1;
s = A*exp(j*2*pi*t*f/Fs);
scomp = zeros(1,2*N);
scomp(1:2:2*N) = real(s);
scomp(2:2:2*N) = imag(s);
save_raw("twotone.iq16",scomp);
