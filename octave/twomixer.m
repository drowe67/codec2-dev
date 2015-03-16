% twomixer.m
% David Rowe Feb 2015
%
% test program for combining signals from two mixers

fc=1024;
Fs=8192;
[b a] = cheby2(6,40,[200 3000]/(Fs/2));

w1 = 2*pi*(fc/2)/(Fs/2);
w2 = 2*pi*(3*fc/2)/(Fs/2);
w3 = 2*pi*(5*fc/2)/(Fs/2);
wc = 2*pi*fc/(Fs/2);

t = 0:(Fs-1);
s1 = cos(w1.*t);
s2 = cos(w2.*t);
s3 = cos(w3.*t);
s = s1 + s2 + s3;

lo1 = cos(wc*t);
lo2 = cos(2*wc*t);
out1 = filter(b,a,s .* lo1);
out2 = filter(b,a,s .* lo2);

S = fft(s);

figure(1)
subplot(211)
plot(abs(S));
subplot(212)
plot(angle(S))

figure(2)
subplot(211)
plot(abs(fft(out1)))
subplot(212)
plot(abs(fft(out2)))
