% nf_from_gr.m
% David Rowe Mar 2016

% TODO: [ ] synthetic signal to test
%       [ ] if sine wave power goes up what happens?

p = load_comp("~/Desktop/neg100dBm_2MHz.bin");
pn = load_comp("~/Desktop/nosignal_2MHz.bin");

Pin = -100;             % rx power input of p in dBm 
Fs = 2E6;               % sample rate in Hz
st = 100E3; en = 300E3; % frequencies window to use
 
P = fft(p100(1:Fs));
N = fft(pns(1:Fs));

PdB = 10*log10(abs(P));
NdB = 10*log10(abs(N));

figure(1)
clf;

subplot(211)
plot(st:en, PdB(st:en));

subplot(212)
plot(st:en, NdB(st:en));

P = 10*log10(P_pwr);
No = 10*log10(No_pwr);
G = P - Pin;
NF  = No - G + 174;
printf("Pin: %4.1f  Pout: %4.1f  G: %4.1f  NF: %3.1f\n", Pin, P, G, NF);
