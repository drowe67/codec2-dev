% phase_nosie.m
% David Nov 2019
%
% Close-in look at phase noise, feed in a off-air sample file of a sine wave

% TODO
%   [ ] estimate freq offset automatically
%   [ ] labels on plots
%   [ ] write float file we can use to run channel simulator
%   [ ] can we estimate bandwith, or some other statistic?  Variance?

function phase_noise(file_name)
  Fs = 8000;
  s = load_raw(file_name);
  S = abs(fft(s(1:Fs).*hanning(Fs)));
  [mx mx_bin] = max(S);
  ftone = mx_bin-1
  
  figure(1); clf;
  plot(20*log10(S(1:Fs/2)))
  title('Input Spectrum');
  
  % downshift to baseband and LPF
  sbb = s' .* exp(-j*(1:length(s))*2*pi*ftone/Fs);
  sbb_lpf = filter(fir1(100,0.1),1,sbb);

  % estimate and remove fine freq offset
  st = Fs; en = 4*Fs;
  phase = unwrap(angle(sbb_lpf(st:en)));
  fine_freq = mean(phase(2:end) - phase(1:end-1));
  sbb_lpf_fine = sbb_lpf .* exp(-j*(1:length(sbb_lpf))*fine_freq);

  figure(2); clf;
  subplot(211)
  plot(real(sbb_lpf(st:en)));
  title('real and imag')
  subplot(212)
  plot(imag(sbb_lpf_fine(st:en)));
  
  figure(3); clf;
  plot(sbb_lpf_fine(st:en))
  title('Polar phase trajectory');

  figure(4); clf;
  S2 = fftshift(fft(sbb_lpf_fine(Fs:Fs*11)));
  middle = round(length(S2)/2)
  [mx mx_bin] = max(abs(S2))
  S2dB = 20*log10(abs(S2));
  mxdB = 10*ceil(max(S2dB)/10);
  x = -10:0.1:10;
  plot(x,S2dB(mx_bin-100:mx_bin+100));
  axis([-10 10 mxdB-40 mxdB])
  title('Close in Phase Spectrum');
  xlabel('Freq (Hz)');
  grid;
  
  figure(5); clf;
  plot(unwrap(angle(sbb_lpf_fine(st:en))))
  title('Unwrapped Phase');
  xlabel('Time (samples)')
  ylabel('Phase (radians)')
  
  figure(6); clf;
  rate_of_change_Hz = (phase(2:end) - phase(1:end-1))*Fs/(2*pi);
  plot(medfilt1(rate_of_change_Hz,1000))
  title('Instantaneous rate of change (Hz)');
  xlabel('Time (samples)')
  ylabel('Freq (Hz)')
end
