% bandpasssampling.m
% David Rowe 23 Feb 2015
%
% Band pass sampling example

graphics_toolkit ("gnuplot");

t = 0:1E-3:1-1E-3;
f1 = 5;
f2 = 105;

x = 1:10:length(s1);

s1 = cos(2*pi*f1*t);
s1_sampled = s1(x);

s2 = cos(2*pi*f2*t);
s2_sampled = s2(x);

figure(1)
subplot(211)
plot(t,s1)
title('5Hz signal sampled at 100 Hz');
subplot(212)
stem(x*1E-3, s1_sampled,'r')
xlabel('Time (s)');

figure(2)
subplot(211)
plot(t,s2)
title('105Hz signal sampled at 100 Hz');
subplot(212)
stem(x*1E-3, s2_sampled,'r')
xlabel('Time (s)');
