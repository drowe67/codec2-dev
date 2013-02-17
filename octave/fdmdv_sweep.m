% fdmdv_sweep.m
% David Rowe Feb 2013
% Produces a raw file that sweeps between 1000 and 2000 Hz to test freq
% response of transmitters.

secs=10;
fmin=1000;
fmax=2000;
Fs=8000;
rms = 4200;  % roughly RMS value of fdmdv signal
amp = sqrt(2)*rms;
nsamples=Fs*secs;
theta = 0;
s=zeros(1,nsamples);

for i=1:nsamples
  f(i) = fmin + i*(fmax-fmin)/nsamples;
  w = 2*pi*f(i)/Fs;  
  theta += w;
  theta -= 2*pi*floor(theta/(2*pi));
  s(i) = amp*cos(theta);  
end

figure(1)
clf
plot(s(1:100));
fout = fopen("1k_2k_sweep.raw", "wb");
fwrite(fout, s, "short");
fclose(fout);

