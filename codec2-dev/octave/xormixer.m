% xormixer.m
% David Rowe Sep 2015
%
% Testing xor gate as a mixer for constant amplitude
% modulation schemes

n = 1024;
carrier = modulation = zeros(1,n);

Tc = 4; % carrier period
for i=1:Tc:n
  carrier(i:i+Tc/2-1) = 1;
end

Tm = 32; % modulation signal period
for i=1:Tm:n
  modulation(i:(i+Tm/2-1)) = 1;
end

%carrier = carrier .* hanning(n)';
%modulation = modulation .* hanning(n)';
mixer = xor(carrier,modulation) .* hanning(n)';

figure(1);
clf
subplot(311)
plot(abs(fft(carrier)))
axis([1 n 0 n/2]);

subplot(312)
plot(abs(fft(modulation)))
axis([1 n 0 n/2]);

subplot(313)
plot(abs(fft(mixer)))
axis([1 n 0 n/2]);

