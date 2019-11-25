% phase_nosie.m
% David Nov 2019
%
% Close-in look at phase noise, feed in a off-air sample file of a sine wave

function phase_noise(file_name)
  Fs = 8000;
  s = load_raw(file_name);
  S = abs(fft(s(1:Fs).*hanning(Fs)));
  [mx mx_bin] = max(S);
  ftone = mx_bin-3
  
  % downshift to baseband and LPF
  sbb = s' .* exp(-j*(1:length(s))*2*pi*ftone/Fs);
  sbb_lpf = filter(fir1(100,0.1),1,sbb);

  figure(1); clf;
  plot(20*log10(S(1:Fs/2)))
  figure(2); clf;
  subplot(211)
  st = Fs; en = 4*Fs;
  plot(real(sbb_lpf(st:en)));
  subplot(212)
  plot(imag(sbb_lpf(st:en)));
  figure(3); clf;
  plot(sbb_lpf(st:en))
  figure(4); clf;
  S2 = fftshift(fft(sbb_lpf(Fs:Fs*11)));
  middle = length(S2)/2
  [mx mx_bin] = max(abs(S2))
  S2dB = 20*log10(abs(S2));
  mxdB = 10*ceil(max(S2dB)/10);
  x = -10:0.1:10;
  plot(x,S2dB(middle-100:middle+100));
  axis([-10 10 mxdB-40 mxdB])
  figure(5); clf;
  plot(unwrap(angle(sbb_lpf(st:en))))
end
