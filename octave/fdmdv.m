% fdmdv.m
%
% Functions that implement a Frequency Divison Multiplexed Modem for
% Digital Voice (FMDV)over HF channels.
%
% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

% reqd to mak sure we get same random bits at mod and demod

rand('state',1); 
randn('state',1);

% Constants

global Fs = 8000;      % sample rate in Hz
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

function tx_symbols = bits_to_qpsk(prev_tx_symbols, tx_bits, modulation)
  global Nc;
  global Nb;

  % re-arrange as (Nc,Nb) matrix

  tx_bits_matrix = zeros(Nc,Nb);
  tx_bits_matrix(1:Nc,1) = tx_bits(1:Nb:Nb*Nc);
  tx_bits_matrix(1:Nc,2) = tx_bits(2:Nb:Nb*Nc);

  if (strcmp(modulation,'dqpsk')) 
    % map to (Nc,1) DQPSK symbols

    for c=1:Nc
      msb = tx_bits_matrix(c,1); lsb = tx_bits_matrix(c,2);

      if ((msb == 0) && (lsb == 0))
	  tx_symbols(c) = prev_tx_symbols(c);
      endif  
      if ((msb == 0) && (lsb == 1))
         tx_symbols(c) = j*prev_tx_symbols(c);
      endif  
      if ((msb == 1) && (lsb == 0))
         tx_symbols(c) = -prev_tx_symbols(c);
      endif  
      if ((msb == 1) && (lsb == 1))
         tx_symbols(c) = -j*prev_tx_symbols(c);
      endif 
  end
  else
    % QPSK mapping
    tx_symbols = -1 + 2*tx_bits_matrix(:,1) - j + 2j*tx_bits_matrix(:,2);
  endif

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

  % add centre pilot tone at twice amplitude of other tones

  c = Nc+1;
  for i=1:M
    phase_tx(c) = phase_tx(c) * freq(c);
    tx_fdm(i) = tx_fdm(i) + 1/sqrt(2)*phase_tx(c);
  end
 
  % Scale such that total Carrier power C of real(tx_fdm) = Nc
  % we return the complex (single sided) signal to make frequency
  % shifting for the purpose of testing easier

  tx_fdm = sqrt(2)*tx_fdm;
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


% Estimate frequency offset of FDM signal by peak picking pilot
% tone at Fcentre

function foff = rx_est_freq_offset(rx_fdm)
  global M;
  global Nc;
  global Fs;
  global freq;
  global phase_rx;
  global Npilotbaseband;
  global Npilotlpf;
  global Npilotcoeff;
  global pilot_baseband;
  global pilot_lpf;
  global pilot_coeff;

  % down convert latest M samples of pilot

  pilot_baseband(1:Npilotbaseband-M) = pilot_baseband(M+1:Npilotbaseband);
  c = Nc+1;
  for i=1:M
    phase_rx(c) = phase_rx(c) * freq(c);
    pilot_baseband(Npilotbaseband-M+i) = rx_fdm(i) * phase_rx(c)';
  end

  % LPF cutoff 200Hz

  pilot_lpf(1:Npilotlpf-M) = pilot_lpf(M+1:Npilotlpf);
  j = 1;
  for i = Npilotlpf-M+1:Npilotlpf
    pilot_lpf(i) = pilot_baseband(j:j+Npilotcoeff) * pilot_coeff';
    j++;
  end

  % decimate to improve DFT resolution, window and DFT

  Mpilot = Fs/(2*200);  % calc decimation rate given new sample rate is twice LPF freq
  Mpilotfft = 256;
  s = pilot_lpf(1:Mpilot:Npilotlpf) .* hanning(Npilotlpf/Mpilot)';
  S = abs(fft(s, Mpilotfft));

  % peak pick and convert to Hz

  [x ix] = max(S);
  r = 2*200/Mpilotfft;     % maps FFT bin to frequency in Hz

  if ix > Mpilotfft/2
    foff = (ix - Mpilotfft - 1)*r;
  else
    foff = (ix - 1)*r;
  endif

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

  rx_timing = angle(x)*M/(2*pi) + M/4;
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
% Maybe phase is coherent over a couple of symbols in HF channel,not
% sure but it's worth 3dB so worth experimenting or using coherent as
% an option.

function rx_phase = rx_est_phase(prev_rx_symbols, rx_symbols)

  % modulation strip

  rx_phase = angle(sum(rx_symbols .^ 4))/4;
 
endfunction


% convert symbols back to an array of bits

function rx_bits = qpsk_to_bits(prev_rx_symbols, rx_symbols, modulation)
  global Nc;
  global Nb;
  global Nb;

  if (strcmp(modulation,'dqpsk')) 
    % extra 45 degree clockwise lets us use real and imag axis as
    % decision boundaries

    phase_difference = rx_symbols .* conj(prev_rx_symbols) * exp(j*pi/4);
  
    % map (Nc,1) DQPSK symbols back into an (1,Nc*Nb) array of bits

    for c=1:Nc
      d = phase_difference(c);
      if ((real(d) >= 0) && (imag(d) >= 0))
         msb = 0; lsb = 0;
      endif  
      if ((real(d) < 0) && (imag(d) >= 0))
         msb = 0; lsb = 1;
      endif  
      if ((real(d) < 0) && (imag(d) < 0))
         msb = 1; lsb = 0;
      endif
      if ((real(d) >= 0) && (imag(d) < 0))
         msb = 1; lsb = 1;
      endif
      rx_bits(2*(c-1)+1) = msb;
      rx_bits(2*(c-1)+2) = lsb;
    end
  else
    % map (Nc,1) QPSK symbols back into an (1,Nc*Nb) array of bits

    rx_bits(1:Nb:Nc*Nb) = real(rx_symbols) > 0;
    rx_bits(2:Nb:Nc*Nb) = imag(rx_symbols) > 0;
  endif

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

  ber = bit_errors/Ntest_bits;
  
  sync = 0;
  if (ber < 0.1)
    sync = 1;
  endif
endfunction


% Initialise ----------------------------------------------------

global tx_filter_memory = zeros(Nc, Nfilter);
global rx_filter_memory = zeros(Nc, Nfilter);

% phasors used for up and down converters

global freq = zeros(Nc+1,1);;
for c=1:Nc/2
  carrier_freq = (-Nc/2 - 1 + c)*Fsep + Fcentre;
  freq(c) = exp(j*2*pi*carrier_freq/Fs);
end
for c=Nc/2+1:Nc
  carrier_freq = (-Nc/2 + c)*Fsep + Fcentre;
  freq(c) = exp(j*2*pi*carrier_freq/Fs);
end

freq(Nc+1) = exp(j*2*pi*Fcentre/Fs);

global phase_tx = ones(Nc+1,1);
global phase_rx = ones(Nc+1,1);

% Freq offset estimator states

global Npilotcoeff    = 30;                             % number of pilot LPF coeffs
global pilot_coeff    = fir1(Npilotcoeff, 200/(Fs/2))'; % 200Hz LPF
global Npilotbaseband = Npilotcoeff + M;                % number of pilot baseband samples reqd for pilot LPF
global Npilotlpf      = 4*M;                            % number of samples we DFT pilot over, pilot est window
global pilot_baseband = zeros(1, Npilotbaseband);       % pilot baseband samples
global pilot_lpf      = zeros(1, Npilotlpf);            % LPC pilot samples

% Timing estimator states

global rx_filter_mem_timing = zeros(Nc, Nt*P);
global rx_baseband_mem_timing = zeros(Nc, Nfiltertiming);

% Test bit stream state variables

global Ntest_bits = Nc*Nb*4;     % length of test sequence
global current_test_bit = 1; 
global test_bits = rand(1,Ntest_bits) > 0.5;
global rx_test_bits_mem = zeros(1,Ntest_bits);

