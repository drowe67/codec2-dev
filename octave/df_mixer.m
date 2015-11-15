% df_mixer.m
%
% David Rowe October 2015
%
% Experimental direction finding using a mixer to FDM signals from two
% antennas.  This simulation tests estimation of phase and freq of 
% mixer LO.
%
% [ ] any block size, dft size, Fs
% [ ] stats for repeated tests
% [ ] add modulation

fm;

function tx = fm_mod
  fm_states.Fs = 192E3;  
  fm_states.fm_max = 3E3;
  fm_states.fd = 5E3;
  fm_states.fc = 48E3;
  fm_states.pre_emp = 0;
  fm_states.de_emp  = 0;
  fm_states.Ts = 1;
  fm_states.output_filter = 1;
  fm_states = analog_fm_init(fm_states);

  nsam = fm_states.Fs;
  t = 0:(nsam-1);
  fm = 1000; wm = 2*pi*fm/fm_states.Fs;
  mod = sin(wm*t);
  tx = analog_fm_mod(fm_states, mod);
endfunction

randn('state',1);
more off;

f     = 48E3;          % signal frequency
flo   = 32.768E3;      % LO frequency
m     = 10^(-12/20);   % modulation index, how far each sideband is down wrt carrier (0.1 == 20dB)

phi   = pi/4;          % phase difference between antennas

theta = 0;             % phase term constant to ant1 and ant2
alpha = 0;             % LO phase

Fs    = 192E3;         % sample rate
CNodB = 40;            % C/No of carrier, which dominates power
                       % Around C/No = 50dB is min for FM (12dB-ish SINAD)
N     = Fs;            % processing block size

% BW is Fs Hz.  No=(total noise power Nt)/Fs. N=noise power=variance of noise source
% C = carrier power = 1
% CNo = C/No = 1/(Nt/Fs), therefore var = Nt = Fs/CNo

CNo = 10^(CNodB/10);
var = Fs/CNo;

t=0:N-1;

w   = 2*pi*f/Fs;
wlo = 2*pi*flo/Fs;

ntrials = 5;
phi_est = zeros(1,ntrials);

for nn=1:ntrials
  theta = 2*pi*rand(1,1);
  alpha = 2*pi*rand(1,1);

  % tx signal

  %tx = exp(j*w*t); 
  tx = fm_mod; 

  % simulate rx signals at antenna 1 and 2 by adding noise.  Note signal
  % at antenna 2 is phase shifted phi due to path difference. Both signals
  % experience a time of flight path different theta.

  n = sqrt(var/2)*randn(1,N) + j*sqrt(var/2)*randn(1,N);
  ant1 = tx*exp(j*(theta+phi)) + n;
  ant2 = tx*exp(j*(theta)) + n;

  % ant2 passed through mixer with local osc frequency wlo and phase alpha

  ant2_mix = 2.0*m*(ant2 .* cos(t*wlo + alpha));

  % The SDR recieves the sum of the two signals

  rx = ant1 + ant2_mix;

  % Lets go tow ork on it and see if we can recover phi, the phase
  % difference between ant1 and ant2

  Rx = (1/N)*fft(rx,N);

  % The peak is the carrier location

  [sam1_mag sam1_ind] = max(abs(Rx));

  sam1_freq = sam1_ind*Fs/N;

  % lets find and sample the peaks either side from the "sidebands"

  lo_offset_ind = flo*N/Fs;
  st = round(0.9*(sam1_ind + lo_offset_ind));
  en = round(1.1*(sam1_ind + lo_offset_ind));
  [sam2_mag sam2_ind] = max(abs(Rx(st:en)));
  sam2_ind += st - 1;
  sam2_freq = sam2_ind*Fs/N;

  st = round(0.9*(sam1_ind - lo_offset_ind));
  en = round(1.1*(sam1_ind - lo_offset_ind));
  [sam3 sam3_ind] = max(abs(Rx(st:en)));
  sam3_ind += st - 1;
  sam3_freq = sam3_ind*Fs/N;

  % Now use the math to find the unknowns ....

  sam1 = Rx(sam1_ind);
  sam2 = Rx(sam2_ind);
  sam3 = Rx(sam3_ind);

  phi_est(nn) = angle(sam1) - (angle(sam2) + angle(sam3))/2;

  printf("carrier: %6.0f Hz LSB: %6.0f Hz USB: %6.0f Hz phi_est: %4.3f\n", sam1_freq, sam3_freq, sam2_freq, phi_est(nn));

  t += N;
end

% some plots ....

figure(1);
clf;
plot((1:N)*Fs/N, 20*log10(abs(Rx)));
axis([1 Fs -60 0])
title('Rx signal at SDR input');

figure(2)
clf
plot(exp(j*phi_est),'+');
axis([-1.5 1.5 -1.5 1.5])

