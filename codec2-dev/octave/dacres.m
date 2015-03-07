% dacres.m
% David Rowe 18 Feb 2015
% Brady O'Brien 5 Mar 2015
% DAC upconversion simulation

graphics_toolkit ("gnuplot");


M = 25;				%interpolation
fs = 2E6/M;			%Input sampling frequency
fb = 7E5;			%Bandpass frequency

f1 = fs/4;
f2 = f1 + 4E3;
f3 = f1 - 4E3;
t = (0:(fs-1));

beta1 = 0.999;
b1x = -2*sqrt(beta1)*cos(2*pi*fb/(fs*M))
beta2 = 1 - (1-0.999)*M;

s1 = [fs zeros(1,fs-1)];       % noise floor, continuous interferers 
s2 = 100*4*cos(t*2*pi*f2/fs);  % wanted signal 40dB above interferers
s3 = 100*3*cos(t*2*pi*f3/fs);  % another wanted frequency
s4 = 100*4*cos(t*2*pi*f1/fs);  % wanted signal at center frequency
s = s1 + s2 + s3 + s4;


s2 = filter([1 0 beta2],1,s);  %pre interpolation notch filter to equalize bandpass after interpolation
s3 = zeros(1,length(s2)*M);    %interpolate by zero-stuffing
s3(1:M:length(s2)*M) = s2;
s4 = filter(1,[1 b1x beta1],s3); %select wanted signal

figure(1)
subplot(211)
plot(20*log10(abs(fft(s)/fs)))
grid
axis([0 fs/2 -10 50])
title('Output from modem');
subplot(212)
plot(20*log10(abs(fft(s2/fs))))
grid
axis([0 fs/2 -10 70])
title('After Pre-eq');

figure(2)
subplot(211)
plot(20*log10(abs(fft(s3)/fs)))
grid
axis([0 (fs/2)*M -10 80])
title('After interpolation');
subplot(212)
plot(20*log10(abs(fft(s4)/fs)))
grid
title('After bandpass');
axis([0 (fs/2)*M -10 80])
