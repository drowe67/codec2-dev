% test_ftrack.m
%
% David Rowe May 2015
%
% Octave script that to test and develop symbol rate frequency offset
% tracking
%
% [X] add Es/No noise and measure BER
% [X] QPSK mod of random symbols
% [X] tracking 1 Hz/s freq drift shiel staying in lock
% [X] simulate R/N filter delay through filter

rand('state',1); 
randn('state',1);
graphics_toolkit ("gnuplot");

Fs = 8000;        % sample rate
Rs = 50;          % QPSK modem symbol rate
Nbits = 5000;     % number of bits to simulute over
Nsymb = Nbits/2;  
foff = 1;         % initial freq offset
dfoff = -1/Fs;    % rate of change of frequency offset in Hz/sample
M = Fs/Rs;        % oversampling rate

EsNodB = 8;       % noise level
EsNo = 10^(EsNodB/10);    
variance = 1/EsNo;

Nsym = 5;         % simulated delay of root nyquist filter in symbols

% control system constants, these can be tweaked

beta = 0.005;
g    = 0.2;

% init state variables

phase_ch = 1;
prev_rx_symb = 1;
ftrack = 0;
filt = 0;

%  QPSK symbols ------> x -------> x ----------G(z) ------------\
%                       |          |                            |
%               exp(j*w_off*n)     |                            |
%                                  |                            |
%                         exp(-j*w_track*n)                     |
%                                  \____________________________/
%
% for small phase/freq errors, where phase/freq detector is linear:
%
%              -1           beta               1            1
%  G(z) = 1 - z     g  ---------------      --------     -------
%                                     -1          -1           -1
%                      1 - (1-beta)z         1 - z        1 - z 
%
%          (1)              (2)               (3)          (4)
%
% (1) - differentiator that gives us an error based on frequency difference
% (2) - first order IIR loop filter  
% (3) - integrator I added so loop could track ramp inputs (i.e. freq drift)
% (4) - frequency to phase integrator implicit in VCO
%
% Note (1) and (3) cancel out, H(z) = 1, so we get a 2nd order loop.  We need to keep (1)
% in the time domain implementation to make the modulation stripper work, it likes angles
% near 0.
%        
% Transfer function, derived with pencil and paper:
% 
%     G(z)         g*beta
% ------------ = ----------  * ------------------------------------
% 1 + G(z)H(z)   1 + g*beta         2 - beta   -1    1 - beta    -2
%                              1 - ---------- z   + ----------- z
%                                  1 + g*beta        1 + g*beta

% Second order IIR filter coeffs

a1 = (2-beta)/(1+g*beta);
a2 = (1-beta)/(1+g*beta);

% Which can be used to work out the polar coordinates of the pole

gamma = sqrt(a2);           % radius (distance from origin)
w = acos(a1/(2*gamma));     % angular frequency

printf("2nd order loop w: %f   gamma: %f\n", w, gamma);


function symbol = qpsk_mod(two_bits)
    two_bits_decimal = sum(two_bits .* [2 1]); 
    switch(two_bits_decimal)
        case (0) symbol =  1;
        case (1) symbol =  j;
        case (2) symbol = -j;
        case (3) symbol = -1;
    endswitch
endfunction


% modulator

tx_bits = rand(1, Nbits) > 0.5;
tx_symb = zeros(1, Nsymb);
for i=1:Nsymb
  tx_symb(i) = qpsk_mod(tx_bits(2*(i-1)+1:2*i));
end

% run rest of simulation symbol by symbol 

rx_symb = zeros(1, Nsymb);
rx_filt = zeros(1,Nsym);

diff_log = [];
foff_log = [];
ftrack_log = [];
mod_strip_log = [];

for i=1:Nsymb

  % channel

  foff_rect = exp(j*2*pi*(foff-ftrack)*M/Fs);
  foff += M*dfoff;
  phase_ch *= foff_rect;
  rx_symb(i) = tx_symb(i) * phase_ch;
  noise = sqrt(variance*0.5)*(randn(1,1) + j*randn(1, 1));
  rx_symb(i) += noise;

  % simulate delay of root nyquist filter

  rx_filt(1:Nsym-1) = rx_filt(2:Nsym);
  rx_filt(Nsym) = rx_symb(i);

  % freq tracking loop

  diff = rx_filt(1) .* conj(prev_rx_symb);
  prev_rx_symb = rx_filt(1);

  % 4th power strips QPSK modulation, by multiplying phase by 4
  % Using the abs value of the real coord was found to help 
  % non-linear issues when noise power was large

  mod_strip = diff.^4;
  mod_strip = abs(real(mod_strip)) + j*imag(mod_strip);

  % loop filter made up of 1st order IIR plus integrator plus gain
  % term.  Integrator was found to be reqd to track ramp inputs
  % (i.e. freq drift)

  filt = ((1-beta)*filt + beta*angle(mod_strip));
  ftrack += g*filt;

  diff_log = [diff_log diff];
  foff_log = [foff_log foff];
  ftrack_log = [ftrack_log ftrack];
  mod_strip_log = [mod_strip_log mod_strip];
end

% plot results

figure(1)
clf
plot(angle(mod_strip_log))
title('mod stripping phase');

figure(2);
clf
plot(mod_strip_log,'+')
axis([-2 2 -2 2])
title('mod stripping scatter');

figure(3)
clf
subplot(211)
plot(real(mod_strip_log))
title('mod stripping real and imag');
subplot(212)
plot(imag(mod_strip_log))

figure(4)
clf
subplot(211)
plot(foff_log,';foff;')
hold on
plot(ftrack_log,'g+;ftrack;')
hold off;
legend("boxoff");  
ylabel('Freq Offset Hz');

subplot(212)
plot(foff_log-ftrack_log,';foff-ftrack;' )
ylabel('Freq Offset Tracking Error Hz');
xlabel('symbols')
legend("boxoff");  

figure(5)
freqz(g*beta/(1+g*beta), [1 -2*gamma*cos(w) gamma*gamma])
