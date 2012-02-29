% fdmdv.m
%
% Octave implementation of a Frequency Divison Multiplexed Modem for
% Digital Voice (FMDV)over HF channels.
%
% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

clear all;
rand('state',1); 
randn('state',1);
 
global Fs = 2000;      % sample rate in Hz
global T  = 1/Fs;      % sample period in seconds
global Rs = 50;        % symbol rate in Hz
global Ts = 1/Rs;      % symbol period in seconds
global Nc = 14;        % number of carriers
global Nb = 2;         % Bits/symbol for QPSK modulation
global Rb = Nc*Rs*Nb;  % bit rate
global M  = Fs/Rs;     % oversampling factor
global Nsym  = 4;      % number of symbols to filter over
global Fsep  = 75;     % Separation between carriers (Hz)
global Fcentre = 1200; % Centre frequency, Nc/2 below this, N/c above (Hz)
global Nt = 5;         % number of symbols we estimate timing over
global P = 4;          % oversample factor used for rx symbol filtering

% Generate root raised cosine (Root Nyquist) filter ---------------

% thanks http://www.dsplog.com/db-install/wp-content/uploads/2008/05/raised_cosine_filter.m

alpha = 0.5;      % excess bandwidth
n = -Nsym*Ts/2:T:Nsym*Ts/2;
global Nfilter = Nsym*M;
global Nfiltertiming = M+Nfilter+M;

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


% Functions ----------------------------------------------------


% generate Nc QPSK symbols from vector of (1,Nc*Nb) input bits

function tx_symbols = bits_to_qpsk(tx_bits)
  global Nc;
  global Nb;

  % re-arrange as (Nc,Nb) matrix

  tx_bits_matrix = zeros(Nc,Nb);
  tx_bits_matrix(1:Nc,1) = tx_bits(1:Nb:Nb*Nc);
  tx_bits_matrix(1:Nc,2) = tx_bits(2:Nb:Nb*Nc);

  % map to (Nc,1) QPSK symbols

  tx_symbols = -1 + 2*tx_bits_matrix(:,1) - j + 2j*tx_bits_matrix(:,2); 

endfunction


% Given Nc*Nb bits construct M samples (1 symbol) of Nc filtered
% symbols streams

function tx_baseband = tx_filter(tx_symbols)
  global Nc;
  global M;
  global tx_filter_memory;
  global Nfilter;
  global gt_alpha5_root;

  tx_baseband = zeros(Nc,M);

  % tx filter each symbol, generate M filtered output samples for each symbol.
  % Efficient polyphase filter techniques used as tx_filter_memory is sparse

  tx_filter_memory(:,Nfilter) = sqrt(2)/2*tx_symbols;
  for i=1:M
    tx_baseband(:,i) = M*tx_filter_memory(:,M:M:Nfilter) * gt_alpha5_root(M-i+1:M:Nfilter)';
  end
  tx_filter_memory(:,1:Nfilter-M) = tx_filter_memory(:,M+1:Nfilter);
  tx_filter_memory(:,Nfilter-M+1:Nfilter) = zeros(Nc,M);
endfunction


% Construct FDM signal by frequency shifting each filtered symbol
% stream

function tx_fdm = fdm_upconvert(tx_filt)
  global Fs;
  global M;
  global Nc;
  global Fsep;
  global phase_tx;
  global freq;

  tx_fdm = zeros(1,M);

  % Nc/2 tones below centre freq
  
  for c=1:Nc/2
      for i=1:M
        phase_tx(c) = phase_tx(c) * freq(c);
	tx_fdm(i) = tx_fdm(i) + tx_filt(c,i)*phase_tx(c);
      end
  end

  % Nc/2 tones above centre freq  

  for c=Nc/2+1:Nc
      for i=1:M
        phase_tx(c) = phase_tx(c) * freq(c);
	tx_fdm(i) = tx_fdm(i) + tx_filt(c,i)*phase_tx(c);
      end
  end

endfunction


% Frequency shift each modem carrier down to Nc baseband signals

function rx_baseband = fdm_downconvert(rx_fdm)
  global Fs;
  global M;
  global Nc;
  global Fsep;
  global phase_rx;
  global freq;

  rx_baseband = zeros(1,M);

  % Nc/2 tones below centre freq
  
  for c=1:Nc/2
      for i=1:M
        phase_rx(c) = phase_rx(c) * freq(c);
	rx_baseband(c,i) = rx_fdm(i)*phase_rx(c)';
      end
  end

  % Nc/2 tones above centre freq  

  for c=Nc/2+1:Nc
      for i=1:M
        phase_rx(c) = phase_rx(c) * freq(c);
	rx_baseband(c,i) = rx_fdm(i)*phase_rx(c)';
      end
  end
endfunction


% Receive filter each baseband signal at oversample rate P

function rx_filt = rx_filter(rx_baseband)
  global Nc;
  global M;
  global P;
  global rx_filter_memory;
  global Nfilter;
  global gt_alpha5_root;
  global Fsep;

  rx_filt = zeros(Nc,P);

  % rx filter each symbol, generate P filtered output samples for each symbol.
  % Note we keep memory at rate M, it's just the filter output at rate P

  N=M/P;
  j=1;
  for i=1:N:M
    rx_filter_memory(:,Nfilter-N+1:Nfilter) = rx_baseband(:,i:i-1+N);
    rx_filt(:,j) = rx_filter_memory * gt_alpha5_root';
    rx_filter_memory(:,1:Nfilter-N) = rx_filter_memory(:,1+N:Nfilter);
    j+=1;
  end
endfunction


% Estimate optimum timing offset, and symbol receive symbols


function [rx_symbols rx_timing] = rx_est_timing(rx_filt, rx_baseband)
  global M;
  global Nt;
  global rx_filter_mem_timing;
  global rx_baseband_mem_timing;
  global P;
  global Nfilter;
  global Nfiltertiming;
  global gt_alpha5_root;

  % update buffer of Nt rate P filtered symbols

  rx_filter_mem_timing(:,1:(Nt-1)*P) = rx_filter_mem_timing(:,P+1:Nt*P);
  rx_filter_mem_timing(:,(Nt-1)*P+1:Nt*P) = rx_filt;

  % sum envelopes of all carriers

  env = sum(abs(rx_filter_mem_timing(:,:)));
  [n m] = size(env);

  % The envelope has a frequency component at the symbol rate.  The
  % phase of this frequency component indicates the timing.  So work out
  % single DFT at frequency 2*pi/P

  x = env * exp(-j*2*pi*(0:m-1)/P)';

  % map phase to estimated optimum timing instant at rate M
  % the M/2 + 1 part was adjusted by experment, I know not why....

  rx_timing = angle(x)*M/pi + M/2 + 1;
  if (rx_timing > M)
     rx_timing -= M;
  end
  if (rx_timing < -M)
     rx_timing += M;
  end

  % rx_filt_mem_timing contains M + Nfilter + M samples of the
  % baseband signal at rate M this enables us to resample the filtered
  % rx symbol with M sample precision once we have rx_timing

  rx_baseband_mem_timing(:,1:Nfiltertiming-M) = rx_baseband_mem_timing(:,M+1:Nfiltertiming);
  rx_baseband_mem_timing(:,Nfiltertiming-M+1:Nfiltertiming) = rx_baseband;

  % sample right in the middle of the timing estimator window, by filtering
  % at rate M

  s = round(rx_timing) + M;
  rx_symbols = rx_baseband_mem_timing(:,s+1:s+Nfilter) * gt_alpha5_root';

endfunction


% Phase estimation function - probably won't work over a HF channel
% Tries to operate over a sinle symbol but uses phase information from
% all Nc carriers which should increase the SNR of phase estimate.
% Maybe phase is coherent over a coupel of symbols in HF channel,not
% sure but it's worth 3dB so worth experimenting or using coherent as
% an option.

function rx_phase = rx_est_phase(rx_symbols)

  % modulation strip

  rx_phase = angle(sum(rx_symbols .^ 4))/4;
 
endfunction


% convert symbols back to an array of bits

function rx_bits = qpsk_to_bits(rx_symbols)
  global Nc;
  global Nb;

  % map (Nc,1) QPSK symbols back into an (1,Nc*Nb) array of bits

  rx_bits = zeros(1,Nc*Nb);
  rx_bits(1:Nb:Nc*Nb) = real(rx_symbols) > 0;
  rx_bits(2:Nb:Nc*Nb) = imag(rx_symbols) > 0;

endfunction


% returns nbits from a repeating sequence of random data

function bits = get_test_bits(nbits)
  global Ntest_bits;       % length of test sequence
  global current_test_bit; 
  global test_bits;

  for i=1:nbits
    bits(i) = test_bits(current_test_bit++);
    if (current_test_bit > Ntest_bits)
      current_test_bit = 1;
    endif
  end

endfunction


% Accepts nbits from rx and attempts to sync with test_bits sequence.
% if sync OK measures bit errors

function [sync bit_errors] = put_test_bits(rx_bits)
  global Ntest_bits;       % length of test sequence
  global test_bits;
  global rx_test_bits_mem;

  % Append to our memory

  [m n] = size(rx_bits);
  rx_test_bits_mem(1:Ntest_bits-n) = rx_test_bits_mem(n+1:Ntest_bits);
  rx_test_bits_mem(Ntest_bits-n+1:Ntest_bits) = rx_bits;

  % see how many bit errors we get when checked against test sequence

  bit_errors = sum(xor(test_bits,rx_test_bits_mem));

  % if less than a thresh we are aligned and in sync

  ber =  bit_errors/Ntest_bits;
  
  sync = 0;
  if (ber < 0.2)
    sync = 1;
  endif
endfunction


% Initialise ----------------------------------------------------

global tx_filter_memory = zeros(Nc, Nfilter);
global rx_filter_memory = zeros(Nc, Nfilter);

% phasors used for up and down converters

global freq = zeros(Nc,1);;
for c=1:Nc/2
  carrier_freq = (-Nc/2 - 1 + c)*Fsep;
  freq(c) = exp(j*2*pi*carrier_freq/Fs);
end
for c=Nc/2+1:Nc
  carrier_freq = (-Nc/2 - 1 + c)*Fsep;
  freq(c) = exp(j*2*pi*carrier_freq/Fs);
end

global phase_tx = ones(Nc,1);
global phase_rx = ones(Nc,1);

global rx_filter_mem_timing = zeros(Nc, Nt*P);
global rx_baseband_mem_timing = zeros(Nc, Nfiltertiming);

frames = 200;
tx_filt = zeros(Nc,M);
rx_symbols_log = zeros(Nc,1);
rx_phase_log = 0;
rx_timing_log = 0;
tx_pwr = 0;
noise_pwr = 0;
total_bit_errors = 0;
total_bits = 0;

global Ntest_bits = 100;     % length of test sequence
global current_test_bit = 1; 
global test_bits = rand(1,Ntest_bits) > 0.5;
global rx_test_bits_mem = zeros(1,Ntest_bits);

% Eb/No calculations.  We need to work out Eb/No for each FDM carrier.
% Total power is sum of power in all FDM carriers

C = 1;  % power of each FDM carrier (energy/sample)
N = 1;  % total noise power (energy/sample) of noise source before scaling by Ngain

EbNo_dB = 2;

% Eb  = Carrier power * symbol time / (bits/symbol)
%     = C *(Rs/Fs) / 2
Eb_dB = 10*log10(C) + 10*log10(Rs) - 10*log10(Fs);

No_dBHz = Eb_dB - EbNo_dB;

% Noise power = Noise spectral density * bandwidth;
N_dB = No_dBHz + 10*log10(Fs) - 10*log10(2);
Ngain_dB = N_dB - 10*log10(N);
Ngain = 10^(Ngain_dB/20);

% C/No = Carrier Power/noise spectral denity
%      = power per carrier*number of carriers / noise spectral denity
CNo_dB = 10*log10(C)  + 10*log10(Nc) - No_dBHz;

% SNR in equivalent 2400 Hz SSB channel

B = 2400;
SNR = CNo_dB - 10*log10(B);


% Main loop ----------------------------------------------------

for i=1:frames
  tx_bits = get_test_bits(Nc*Nb);
  tx_symbols = bits_to_qpsk(tx_bits);
  tx_baseband = tx_filter(tx_symbols);
  tx_fdm = fdm_upconvert(tx_baseband);
  tx_pwr = 0.9*tx_pwr + 0.1*tx_fdm*tx_fdm'/(M);

  noise = Ngain/sqrt(2)*[randn(1,M) + j*randn(1,M)];
  noise_pwr = 0.9*noise_pwr + 0.1*noise*noise'/M;
  rx_fdm = tx_fdm + noise;

  rx_baseband = fdm_downconvert(rx_fdm);
  rx_filt = rx_filter(rx_baseband);

  [rx_symbols rx_timing] = rx_est_timing(rx_filt, rx_baseband);
  rx_timing_log = [rx_timing_log rx_timing];

  %rx_phase = rx_est_phase(rx_symbols);
  %rx_phase_log = [rx_phase_log rx_phase];
  %rx_symbols = rx_symbols*exp(j*rx_phase);

  rx_symbols_log = [rx_symbols_log rx_symbols];
  rx_bits = qpsk_to_bits(rx_symbols);

  [sync bit_errors] = put_test_bits(rx_bits);
  if (sync == 1)
    total_bit_errors = total_bit_errors + bit_errors;
    total_bits = total_bits + Ntest_bits;
  end

  %if (i > 20)
  %  bit_errors = sum(xor(tx_bits,rx_bits));
  %  total_bit_errors = total_bit_errors + bit_errors;
  %  total_bits = total_bits + Nc*Nb;
  %endif
end

ber = total_bit_errors/total_bits;
printf("Eb/No: %2.2f dB  %d bits  %d errors  Meas BER: %1.4f  Theor BER: %1.4f\n", EbNo_dB, 
      total_bits, total_bit_errors, ber, 0.5*erfc(sqrt(10.^(EbNo_dB/10))));

figure(1)
clf;
[n m] = size(rx_symbols_log);
plot(real(rx_symbols_log(:,20:m)),imag(rx_symbols_log(:,20:m)),'+')
figure(2)
clf;
subplot(211)
plot(rx_timing_log)
subplot(212)
plot(180/pi*unwrap(rx_phase_log))

% TODO
%   + state machine to sync on PR test and measure BER
%   + handling sample slips, extra plus/minus samples
%   + simulating sample clock offsets
%   + timing, frequency offset get function
%   + memory of recent tx and rx signal for spectrogram
%   + scatter diagram get function

% DPSK
% add channel noise, phase offset, frequency offset, timing offset
% fading simulator
% file I/O to test with Codec
% sync recovery time
% dump file type plotting & instrumentation
% determine if error pattern is bursty
% HF channel simulation

% BER issues:
%   QPSK mapping
%   Single sided noise issues
%   interference between carriers due to excess BW
%   Crappy RN coeffs
%   timing recovery off by one
%   Use a true PR sequence
%   Sensitivity to Fs
