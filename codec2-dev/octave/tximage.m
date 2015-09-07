% tximage.m
% David Rpowe Sep 2015

% Say we have two signals from a mixer, at f_lo +/- f_if
% We attenuate one by 20dB with a BPF
% We then pass through a non linear amp
% Note non linearities are tricky to simulate as sign() gives
% us a lot of HF harmonics that aliased back down.
%
% However, it does appear that non-linear PA(sign) only supresses image
% by another 6dB.  So mixing up does have it's problems compared
% to direct generation.  A softer limiter (sqrt) seems to do a better job.

fs = 1E4;
f_rf = 146;
f_if = 10.7;
f_lo = f_rf - f_if;

f_rf1 = f_if + f_lo;
f_rf2 = f_if - f_lo;

a_rf1 = 1;
a_rf2 = 0.1;

t = 1:fs;

sig_rf = a_rf1*cos(2*pi*t*f_rf1/fs) + a_rf2*cos(2*pi*t*f_rf2/fs);
sig_rf_clip = sign(sig_rf).* sig_rf .^ 1/2;
b = fir1(100, f_rf/fs);
sig_rf_clip_filter = filter(b,1,sig_rf_clip);

figure(1);
subplot(211)
plot(sig_rf(1:1000))
subplot(212)
plot(sig_rf_clip_filter(1:1000))
figure(2)
plot(20*log10(abs(fft(sig_rf_clip_filter))))
axis([1 1000 0 80])
grid
