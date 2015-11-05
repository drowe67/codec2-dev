% df_mixer.m
%
% David Rowe October 2015
%
% Experimental direction fing using a mixer to FDM signals from two
% antennas.  This simulation tests estimation of phase and freq of 
% mixer LO.
%
% [ ] any block size, dft size, Fs
% [ ] stats for repeated tests
% [ ] add modulation

randn('state',1);

f     = 48E3;   % signal frequency
flo   = 32E3;   % LO frequency
m     = 0.1;    % modulation index, how far each sideband is down wrt carrier (0.1 == 20dB)
phi   = pi/4;   % LO phase
Fs    = 96E3;   % sample rate
CNodB = 40;     % C/No of carrier, which dominates power
                % Around C/No = 50dB is min for FM (12dB-ish SINAD)
N     = Fs;     % processing block size

% BW is Fs Hz.  No=(total noise power Nt)/Fs. N=noise power=variance of noise source
% C = carrier power = 1
% CNo = C/No = 1/(Nt/Fs), therefore var = Nt = Fs/CNo

CNo = 10^(CNodB/10);
var = Fs/CNo;

t=0:N-1;

% tx signal

tx = exp(j*t*2*pi*f/Fs); 

% simulate rx signals at antenna 1 and 2 by adding noise.  Note signal
% at antenna 2 is phase shifted

ant1 = tx + sqrt(var/2)*randn(1,N) + j*sqrt(var/2)*randn(1,N);
ant2 = tx*exp(j*phi) + sqrt(var/2)*randn(1,N) + j*sqrt(var/2)*randn(1,N);

% ant2 passed through mixer, that attenuates signal quite a bit and
% has poor (0dB) carrier supression at UHF

ant2_mix = m*(ant2 .* (1 + 2*cos(t*2*pi*flo/Fs + phi)));

% sum the two signals, to get the overal signal that the SDR samples

rx = ant1 + ant2_mix;

% Turns out a carrier with two sidebands looks at lot like AM.  That
% was a lucky break! So lets simply envelope detect it.

d = abs(rx);
D = (1/N)*fft(d);

figure(1);
clf;
plot((1:N)*Fs/N, 20*log10(abs((1/N)*fft(rx))));
axis([1 Fs -60 0])
title('Rx signal at SDR input');

figure(2);
clf;
subplot(211)
plot((1:N)*Fs/N, 20*log10(abs(D)))
axis([1 Fs -60 0])
title('AM demodulated');
subplot(212)
plot(d(1:10*Fs/flo))


% search for maxima due to LO/sidebands.  We need a search range to
% avoid big DC component

st = round(0.9*N*flo/Fs);
en = round(1.1*N*flo/Fs);
[mx mx_ind] = max(abs(D(st:en)));
phi_est = angle(D(st-1+mx_ind));
st
en
mx_ind
mx
abs(D(st-1+mx_ind))
D(st-1+mx_ind-1:st-1+mx_ind+1)
printf("N: %d C/No (dB): %2.0f dB C/No (lin): %3.1f var: %f\n", N, CNodB, CNo, var);
printf("phi: %4.3f phi: %4.3f\n", phi, phi_est);


