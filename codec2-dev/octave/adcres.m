% adcres.m
% David Rowe 18 Feb 2015
%
% ADC resamping simulation, IIR tuner development.

% [ ] quantisation of ADC
% [ ] SNR at ADC input, SNR at resampler output
% [X] decimation to 50 kHz
% [X] 40dB ish rejection (test)
% [X] visualise pass band flatness

graphics_toolkit ("gnuplot");

fs = 2E6;
f1 = 500E3;
f2 = f1 + 8E3;
f3 = f1 - 7E3;
f4 = f1 - 207E3;
t = (0:(fs-1));
M = 45;
beta1 = 0.999;
beta2 = 1 - (1-beta1)*M;

s1 = [fs zeros(1,fs-1)];       % noise floor, continuous interferers 
s2 = 100*4*cos(t*2*pi*f2/fs);  % wanted signal 40dB above interferers
s3 = 100*2*cos(t*2*pi*f3/fs);
s4 = 100*2*cos(t*2*pi*f4/fs);  % interferer at same level
s = s1 + s2 + s3 + s4;

s2 = filter(1,[1 0 beta1],s);  % BPF at fs/4
s3 = s2(1:M:length(s2));       % decimate 
s4 = filter([1 0 beta2],1,s3); % flatten filter response again

figure(1)
subplot(211)
plot(20*log10(abs(fft(s)/fs)))
grid
axis([0 fs/2 -10 50])
title('Input to ADC');
subplot(212)
plot(20*log10(abs(fft(s2/fs))))
grid
axis([0 fs/2 -10 70])
title('After BPF');

figure(2)
subplot(211)
plot(20*log10(abs(fft(s3)/fs)))
grid
axis([0 fs/2/M -10 50])
title('After Decimation');
subplot(212)
plot(20*log10(abs(fft(s4)/fs)))
grid
title('After Equaliser');
axis([0 fs/2/M -10 50])
