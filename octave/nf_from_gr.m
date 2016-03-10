% nf_from_gr.m
% David Rowe Mar 2016
%
% Calculate NF from GNU Radio ...IQIQ... (32 bit float) sample files

% Take one sample with a -100dBm input carrier
% Take another sample with no signal (just rx noise)

% HackRF --------------------------

%p = load_comp("~/Desktop/nf/hackrf_100dbm_4MHz.bin");
%pn = load_comp("~/Desktop/nf/hackrf_nosignal_4MHz.bin");
%Fs = 4E6;                 % sample rate in Hz
%st = 180E3; en = 400E3;   % frequencies window to use

% RTL-SDR --------------------------

%p = load_comp("~/Desktop/nf/neg100dBm_2MHz.bin");
%pn = load_comp("~/Desktop/nf/nosignal_2MHz.bin");
%Fs = 2E6;               % sample rate in Hz
%st = 100E3; en = 300E3;   % frequencies window to use

% AirSpy -------------------------

p =  load_comp("~/Desktop/nf/airspy_gain_-100dBm_2.5MHz.bin");
pn = load_comp("~/Desktop/nf/airspy_gain_nosig_2.5MHz.bin");
Fs = 2.5E6;               % sample rate in Hz
st = 795E3; en = 820E3;   % frequencies window to use

% Fun Cube Dongle Pro Plus -------------------------

%p =  load_comp("~/Desktop/nf/fcdpp_100dbm_192khz.bin");
%pn = load_comp("~/Desktop/nf/fcdpp_nosig_192khz.bin");
%Fs = 192E3;               % sample rate in Hz
%st = 25E3; en = 125E3;    % frequencies window to use

Pin = -100;               % rx power input of p in dBm 
 
P = fft(p(1:Fs));
N = fft(pn(1:Fs));

PdB = 10*log10(abs(P));
NdB = 10*log10(abs(N));

figure(1)
clf;

subplot(211)
plot(st:en, PdB(st:en));

subplot(212)
plot(st:en, NdB(st:en));

P_dB = 10*log10(var(P(st:en)));
No_dB = 10*log10(var(N(st:en))/(en-st));
G = P_dB - Pin;
NF  = No_dB - G + 174;
printf("Pin: %4.1f  Pout: %4.1f  G: %4.1f  NF: %3.1f\n", Pin, P_dB, G, NF);
