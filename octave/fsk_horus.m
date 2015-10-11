% fsk_horus.txt
% David Rowe 10 Oct 2015
%
% Experimental near space balloon FSK demodulator
% Assume high SNR, but fades near end of mission can wipe out a few bits
% So low SNR perf not a huge issue
%
% [ ] processing buffers of 1 second
%     + 8000 samples input
%     + keep 30 second sliding window to extract packet from
%     + do fine timing on this
%     [X] estimate frequency of two tones
%         + this way we cope with variable shift and drift
%     [ ] estimate amplitudes and equalise, or limit
% [ ] bit flipping against CRC
% [ ] implement CRC
% [ ] frame sync
% [ ] fine timing
% [ ] compare to fldigi

1;

function states = fsk_horus_init()
  states.Ndft = 8192;
  Fs = states.Fs = 8000;
  N = states.N = Fs;             % processing buffer size, nice big window for f1,f2 estimation
  Rs = states.Rs = 50;
  Ts = states.Ts = Fs/Rs;
  states.nsym = N/Ts;
  Nmem = states.Nmem  = N+2*Ts;    % two symbols memory
  states.f1_dc = zeros(1,Nmem);
  states.f2_dc = zeros(1,Nmem);
  states.P = 4;                  % oversample rate out of filter
  states.nin = N;                % can be N +/- Ts/P = N +/- 40 samples to adjust for sample clock offsets
  states.verbose = 1;
endfunction


% test modulator function

function tx  = fsk_horus_mod(states, tx_bits)
    tx = zeros(states.Ts*length(tx_bits),1);
    tx_phase = 0;
    Ts = states.Ts;
    f1 = 1500; f2 = 1900;

    for i=1:length(tx_bits)
      for k=1:Ts
        if tx_bits(i) == 1
          tx_phase += 2*pi*f1/states.Fs;
        else
          tx_phase += 2*pi*f2/states.Fs;
        end
        tx_phase = tx_phase - floor(tx_phase/(2*pi))*2*pi;
        tx((i-1)*Ts+k) = cos(tx_phase);
      end
    end

endfunction


% down converts and filters FSK audio samples, returning the filter output
% for f1 and f2, ready for fine timing estimation.  The f1 and f2 samples 
% are at a sample rate of 4Rs.
%
% Automagically estimates the frequency of the two tones, or
% looking at it another way, the frequency offset and shift
%
% nin is the number of input samples required by demodulator.  This is
% time varying.  It will nominally be N (8000), and occasionally N +/-
% Ts/P (8040 or 7960).  This is how we compensate for differences between the
% remote tx sample clock and our sample clock.  This function always returns
% N/Ts (50) demodulated bits.  Variable number of input samples, constant number
% of output bits.

function [f1_int f2_int states] = fsk_horus_dc_filter(states, sf)
  N = states.N;
  Ndft = states.Ndft;
  Fs = states.Fs;
  Rs = states.Rs;
  Ts = states.Ts;
  nsym = states.nsym;
  P = states.P;
  nin = states.nin;
  verbose = states.verbose;

  assert(length(sf) == nin);

  % find tone frequency and amplitudes ---------------------------------------------

  h = hanning(nin);
  Sf = fft(sf .* h, Ndft);
  [m1 m1_index] = max(Sf(1:Ndft/2));

  % zero out region around max so we can find second highest peak

  Sf2 = Sf;
  st = m1_index - 100;
  if st < 1
    st = 1;
  end
  en = m1_index + 100;
  if en > Ndft/2
    en = Ndft/2;
  end
  Sf2(st:en) = 0;

  [m2 m2_index] = max(Sf2(1:Ndft/2));

  % f1 always the lower tone

  if m1_index < m2_index
    f1 = (m1_index-1)*Fs/Ndft;
    f2 = (m2_index-1)*Fs/Ndft;
    twist = 20*log10(m1/m2);
  else
    f1 = (m2_index-1)*Fs/Ndft;
    f2 = (m1_index-1)*Fs/Ndft;
    twist = 20*log10(m2/m1);
  end

  states.f1 = f1;
  states.f2 = f2;

  % down convert and filter at 4Rs ------------------------------

  % update filter (integrator) memory by shifing in nin samples
  
  nold = N+2*Ts-nin; % number of old samples we retain

  f1_dc = states.f1_dc; 
  f1_dc(1:nold) = f1_dc(N+2*Ts-nold+1:N+2*Ts);
  f2_dc = states.f2_dc; 
  f2_dc(1:nold) = f2_dc(N+2*Ts-nold+1:N+2*Ts);

  % shift down to around DC

  f1_dc(nold+1:N+2*Ts) = sf' .* exp(-j*(0:nin-1)*2*pi*f1/Fs);
  f2_dc(nold+1:N+2*Ts) = sf' .* exp(-j*(0:nin-1)*2*pi*f2/Fs);

  % save filter (integrator) memory for next time

  states.f1_dc = f1_dc;
  states.f2_dc = f2_dc;

  % integrate over symbol period, which is effectively a LPF, removing
  % the -2Fc frequency image.  Can also be interpreted as an ideal
  % integrate and dump, non-coherent demod.  We run the integrator at
  % PRs = 4Rs (1/4 symbol offsets) to get outputs at a range of different
  % fine timing offsets.  We calculate integrator output over nsym+1 (e.g. 51)
  % symbols so we have extra samples for the fine timing re-sampler at either
  % end of the array.

  rx_bits = zeros(1, (nsym+1)*P);
  for i=1:(nsym+1)*P
    st = 1 + (i-1)*Ts/P;
    en = st+Ts-1;
    f1_int(i) = sum(f1_dc(st:en));
    f2_int(i) = sum(f2_dc(st:en));
  end

  % fine timing estimation -----------------------------------------------

  % Non linearity has a spectral line at Rs, with a phase
  % related to the fine timing offset.  See:
  %   http://www.rowetel.com/blog/?p=3573 
  % We have sampled the integrator output at Fs=P samples/symbol, so
  % lets do a single point DFT at w = 2*pi*f/Fs = 2*pi*Rs/(P*Rs)

  Np = length(f1_int);
  w = 2*pi*(Rs)/(P*Rs);
  x = abs((f1_int-f2_int)) * exp(-j*w*(0:Np-1))';
  norm_rx_timing = angle(x)/(2*pi);
  rx_timing = norm_rx_timing*P;

  % keep within 1..P-1

  if (rx_timing > P-1)
     rx_timing -= P;
  end
  if (rx_timing < 0)
     rx_timing += P;
  end

  % work out how many input samples we need next time, to keep rx_timing inside a 0.5 symbol range

  next_nin = N;
  if rx_timing > 3
     next_nin += Ts/P;
  end
  if rx_timing < 1;
     next_nin -= Ts/P;
  end
  states.nin = next_nin;

  % Re sample integrator outputs using fine timing estimate

  low_sample = floor(rx_timing);
  fract = rx_timing - low_sample;
  high_sample = ceil(rx_timing);

  if verbose
    printf("rx_timing: %3.2f low_sample: %d high_sample: %d fract: %3.3f nin_next: %d\n", rx_timing, low_sample, high_sample, fract, next_nin);
  end

  f1_int_resample = zeros(1,nsym);
  f2_int_resample = zeros(1,nsym);
  for i=1:nsym
    st = (i-1)*P + 1;
    f1_int_resample(i) = f1_int(st+low_sample)*(1-fract) + f1_int(st+high_sample)*fract;
    f2_int_resample(i) = f2_int(st+low_sample)*(1-fract) + f2_int(st+high_sample)*fract;
  end

  states.f1_int_resample = f1_int_resample;
  states.f2_int_resample = f2_int_resample;
endfunction


% demo script --------------------------------------------------------

more off
states = fsk_horus_init();
N = states.N;
P = states.P;
Rs = states.Rs;

EbNodB = 100;
EbNo = 10^(EbNodB/10);
variance = states.Fs/(states.Rs*EbNo);

%s = load_raw("~/Desktop/vk5arg-3.wav");

tx_bits = round(rand(1, N*10/states.Ts));
%tx_bits = zeros(1,3*N/states.Ts);
%tx_bits(1:2:length(tx_bits)) = 1;
s  = fsk_horus_mod(states, tx_bits);
s  += sqrt(variance/2)*randn(length(s),1);

timing_offset = round(0.3*states.Ts);
st = 1 + timing_offset;
for f=1:8

  % extract nin samples from input stream

  nin = states.nin;
  en = st + states.nin - 1;
  sf = s(st:en);
  st += nin;

  % demodulate to stream of bits
  [f1_int f2_int states] = fsk_horus_dc_filter(states, sf);
end

figure(1)
plot(abs(f1_int),'+')
hold on;
plot(abs(f2_int),'g+')
hold off;

figure(2);
plot(abs(states.f1_int_resample),'+')
hold on;
plot(abs(states.f2_int_resample),'g+')
hold off;

