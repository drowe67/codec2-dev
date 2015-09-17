% adc_sfdr_ut.m
% David Rowe Aug 2015
%
% Processes data collected from STM32F4 or SFDR testing of ADC

s = load_raw("~/stlink/adc.raw");
Fs = 2E6;
N = 1024;
num_frames = length(s)/N;
h = hanning(N);
XdB = zeros(N/2,1);

for i=1:num_frames
  x = s((i-1)*N+1:i*N);
  X = fft(x .* h);
  XdB += 20*log10(abs(X(1:N/2)));
end

XdB /= num_frames;
XdB -= max(20*log10(N));

figure(1)
clf
plot((0:N/2-1)*Fs/(1000*N), XdB)
grid
ylabel('Amplitude dB')
xlabel('Frequency (kHz)');
axis([0 Fs/(2*1000) -30 80])

figure(2)
clf
plot(s)


