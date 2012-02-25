% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Octave implementation of a G3PLX style FDMDV HF modem.

clear all;
rand('state',1); 
randn('state',1);
 
global Fs = 1200;      % sample rate in Hz
global T  = 1/Fs;      % sample period in seconds
global Rs = 50;        % symbol rate in Hz
global Ts = 1/Rs;      % symbol period in seconds
global Nc = 14;        % number of carriers
global Nb = 2;         % Bits/symbol for QPSK modulation
global Rb = Nc*Rs*Nb;  % bit rate
global M  = Fs/Rs;     % oversampling factor
global Nsym  = 8;      % number of symbols to filter over
global Fsep  = 75;     % Separation between carriers (Hz)
global Fcentre = 1200; % Centre frequency, Nc/2 below this, N/c above (Hz)

% Generate root raised cosine (Root Nyquist) filter ---------------
% thanks http://www.dsplog.com/db-install/wp-content/uploads/2008/05/raised_cosine_filter.m

alpha = 0.5;      % excess bandwidth
n = -Nsym*Ts/2:T:Nsym*Ts/2;
global Nfilter = Nsym*M + 1;

sincNum = sin(pi*n/Ts); % numerator of the sinc function
sincDen = (pi*n/Ts);    % denominator of the sinc function
sincDenZero = find(abs(sincDen) < 10^-10);
sincOp = sincNum./sincDen;
sincOp(sincDenZero) = 1; % sin(pix/(pix) =1 for x =0

cosNum = cos(alpha*pi*n/Ts);
cosDen = (1-(2*alpha*n/Ts).^2);
cosDenZero = find(abs(cosDen)<10^-10);
cosOp = cosNum./cosDen;
cosOp(cosDenZero) = pi/4;
gt_alpha5 = sincOp.*cosOp;
Nfft = 4096;
GF_alpha5 = fft(gt_alpha5,Nfft)/M;
% sqrt causes stop band to be amplifed, this hack pushes it down again
for i=1:Nfft
  if (abs(GF_alpha5(i)) < 0.02)
    GF_alpha5(i) *= 0.001;
  endif
end
GF_alpha5_root = sqrt(abs(GF_alpha5)) .* exp(j*angle(GF_alpha5));
ifft_GF_alpha5_root = ifft(GF_alpha5_root);
global gt_alpha5_root = real((ifft_GF_alpha5_root(1:Nfilter)));

% Modulator ----------------------------------------------------

% Given Nc*Nb bits construct M samples (1 symbol) of Nc filtered
% symbols streams

function tx_filt = filter_symbols(tx_bits)
  global Nc;
  global M;
  global tx_filter_memory;
  global Nfilter;
  global gt_alpha5_root;
  global Fsep;

  tx_filt = zeros(Nc,M);

  % generate Nc QPSK symbols from Nc*Nb input bits

  tx_symbols = 1 - 2*tx_bits(:,1) + j - 2j*tx_bits(:,2); 

  % tx filter each symbol, generate M filtered output samples for each symbol

  tx_filter_memory(:,Nfilter) = tx_symbols;
  for i=1:M
    tx_filt(:,i) = M*tx_filter_memory * gt_alpha5_root';
    tx_filter_memory(:,1:Nfilter-1) = tx_filter_memory(:,2:Nfilter);
    tx_filter_memory(:,Nfilter) = zeros(Nc,1);
  end
endfunction

% Construct FDM signal that combines all filtered symbols by
% frequency shifting each filtered symbol stream

function tx_fdm = build_fdm(tx_filt)
  global Fs;
  global M;
  global Nc;
  global Fsep;
  global phase_tx;

  tx_fdm = zeros(1,M);

  % Nc/2 tones below centre freq
  
  for c=1:Nc/2
      carrier_freq = (-Nc/2 - 1 + c)*Fsep;
      for i=1:M
        phase_tx(c) = phase_tx(c) + 2*pi*carrier_freq/Fs;
	tx_fdm(i) = tx_fdm(i) + tx_filt(c,i)*exp(j*phase_tx(c));
      end
  end

  % Nc/2 tones above centre freq  

  for c=Nc/2+1:Nc
      carrier_freq = (-Nc/2 + c)*Fsep;
      for i=1:M
        phase_tx(c) = phase_tx(c) + 2*pi*carrier_freq/Fs;
	tx_fdm(i) = tx_fdm(i) + tx_filt(c,i)*exp(j*phase_tx(c));
      end
  end

endfunction


% Main loop ----------------------------------------------------

global tx_filter_memory = zeros(Nc, Nfilter);

frames = 100;
tx_filt = zeros(Nc,M);
global phase_tx = zeros(Nc,1);
tx_fdm = zeros(1,M);

for i=1:frames
  tx_bits = rand(Nc,Nb) > 0.5; 
  tx_filt = filter_symbols(tx_bits);
  tx_fdm = [tx_fdm build_fdm(tx_filt)];
end

figure(1)
clf;
plot(abs(fft((tx_fdm))))

% test demod

[m n] = size(tx_fdm)
rx_bb = tx_fdm .* exp(-j*(0:n-1)*2*pi*Fsep/Fs);
rx_bb_filt = conv(gt_alpha5_root, rx_bb);

figure(2)
clf;
plot(real(rx_bb_filt))

% add channel noise, phase offset, frequency offset, timing offset

% root nyquist rx filter

% timing recovery - sum envelopes



