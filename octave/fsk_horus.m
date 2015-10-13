% fsk_horus.m
% David Rowe 10 Oct 2015
%
% Experimental near space balloon FSK demodulator
% Assume high SNR, but fades near end of mission can wipe out a few bits
% So low SNR perf not a huge issue
%
% [X] processing buffers of 1 second
%     + 8000 samples input
%     + keep 30 second sliding window to extract packet from
%     + do fine timing on this
%     [X] estimate frequency of two tones
%         + this way we cope with variable shift and drift
%         + starts to lose it at 8 Eb/No = 8db.  Maybe wider window?
%     [X] estimate amplitudes and equalise, or limit
%         + not needed - tones so close in freq unlikely to be significant ampl diff
%           across SSB rx filter
% [X] Eb/No point 8dB, 2% ish
% [X] fine timing and sample slips, +/- 1000ppm (0.1%) clock offset test
% [ ] bit flipping against CRC
% [ ] implement CRC
% [X] frame sync
% [X] compare to fldigi
%     + in AWGN channel 3-4dB improvement.  In my tests fldigi can't decode  
%       with fading model, requires Eb/No > 40dB, this demo useable at Eb/No = 20dB
% [ ] test over range of f1/f2, shifts, timing offsets, clock offsets, Eb/No
%     [X] +/- 1000ppm clock offset OK at Eb/No = 10dB, starts to lose it at 8dB
%     [X] tone freq est starts to lose it at 8dB in awgn.  Maybe wider window?
% [ ] low snr detection of $$$$$$
%     + we might be able to pick up a "ping" at very low SNRs to help find baloon on ground
% [ ] streaming, indicator of audio freq, i.e. speaker output

1;

function states = fsk_horus_init()
  states.Ndft = 8192;
  Fs = states.Fs = 8000;
  N = states.N = Fs;             % processing buffer size, nice big window for f1,f2 estimation
  Rs = states.Rs = 50;
  Ts = states.Ts = Fs/Rs;
  states.nsym = N/Ts;
  Nmem = states.Nmem  = N+2*Ts;    % two symbol memory in down converted signals to allow for timing adj
  states.f1_dc = zeros(1,Nmem);
  states.f2_dc = zeros(1,Nmem);
  states.P = 8;                  % oversample rate out of filter
  states.nin = N;                % can be N +/- Ts/P = N +/- 40 samples to adjust for sample clock offsets
  states.verbose = 0;
  states.phi1 = 0;               % keep down converter osc phase continuous
  states.phi2 = 0;

  % Generate unque word that correlates against the ASCIi "$$$$$" that
  % delimits start and end of frame Note use of zeros in UW as "don't
  % cares", we ignore RS232 start/stop bits.  Not sure this is a good
  % idea, we could include start and stop bits if we like.  Oh Well.

  dollar_bits = fliplr([0 1 0 0 1 0 0]);
  mapped_db = 2*dollar_bits - 1;
  npad = states.npad = 3;   % one start and two stop bits between 7 bit ascii chars
  nfield = states.nfield = 7; % length of ascii character field

  states.uw = [mapped_db zeros(1,npad) mapped_db zeros(1,npad) mapped_db zeros(1,npad)  mapped_db zeros(1,npad) mapped_db zeros(1,npad)];

  states.uw_thresh = 5*7 - 4; % allow 24bit errors when looking for UW

  states.df = 0;
  states.f1 = 0;
  states.f2 = 0;
endfunction


% test modulator function

function tx  = fsk_horus_mod(states, tx_bits)
    tx = zeros(states.Ts*length(tx_bits),1);
    tx_phase = 0;
    Ts = states.Ts;
    Fs = states.Fs;
    f1 = 1000; f2 = 1400;
    df = states.df; % tone freq change in Hz/s

    for i=1:length(tx_bits)
      if tx_bits(i) == 0
        tx_phase_vec = tx_phase + (1:Ts)*2*pi*f1/Fs;
      else
        tx_phase_vec = tx_phase + (1:Ts)*2*pi*f2/Fs;
      end
      tx((i-1)*Ts+1:i*Ts) = 2.0*cos(tx_phase_vec);
      tx_phase -= floor(tx_phase/(2*pi))*2*pi;
      f1 += df*Ts/Fs; f2 += df*Ts/Fs;
    end
endfunction


% Given a buffer of nin input Rs baud FSK samples, returns nsym bits.
%
% Automagically estimates the frequency of the two tones, or
% looking at it another way, the frequency offset and shift
%
% nin is the number of input samples required by demodulator.  This is
% time varying.  It will nominally be N (8000), and occasionally N +/- 
% Ts/2 (e.g. 8080 or 7920).  This is how we compensate for differences between the
% remote tx sample clock and our sample clock.  This function always returns
% N/Ts (50) demodulated bits.  Variable number of input samples, constant number
% of output bits.

function [rx_bits states] = fsk_horus_demod(states, sf)
  N = states.N;
  Ndft = states.Ndft;
  Fs = states.Fs;
  Rs = states.Rs;
  Ts = states.Ts;
  nsym = states.nsym;
  P = states.P;
  nin = states.nin;
  verbose = states.verbose;
  Nmem = states.Nmem;

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

  if verbose
    %printf("f1: %4.0f Hz f2: %4.0f Hz\n", f1, f2);
  end

  % down convert and filter at rate P ------------------------------

  % update filter (integrator) memory by shifting in nin samples
  
  nold = Nmem-nin; % number of old samples we retain

  f1_dc = states.f1_dc; 
  f1_dc(1:nold) = f1_dc(Nmem-nold+1:Nmem);
  f2_dc = states.f2_dc; 
  f2_dc(1:nold) = f2_dc(Nmem-nold+1:Nmem);

  % shift down to around DC, ensuring continuous phase from last frame

  phi1_vec = states.phi1 + (1:nin)*2*pi*f1/Fs;
  phi2_vec = states.phi2 + (1:nin)*2*pi*f2/Fs;

  f1_dc(nold+1:Nmem) = sf' .* exp(-j*phi1_vec);
  f2_dc(nold+1:Nmem) = sf' .* exp(-j*phi2_vec);

  states.phi1  = phi1_vec(nin);
  states.phi1 -= 2*pi*floor(states.phi1/(2*pi));
  states.phi2  = phi2_vec(nin);
  states.phi2 -= 2*pi*floor(states.phi2/(2*pi));

  % save filter (integrator) memory for next time

  states.f1_dc = f1_dc;
  states.f2_dc = f2_dc;

  % integrate over symbol period, which is effectively a LPF, removing
  % the -2Fc frequency image.  Can also be interpreted as an ideal
  % integrate and dump, non-coherent demod.  We run the integrator at
  % rate P (1/P symbol offsets) to get outputs at a range of different
  % fine timing offsets.  We calculate integrator output over nsym+1
  % symbols so we have extra samples for the fine timing re-sampler at either
  % end of the array.

  rx_bits = zeros(1, (nsym+1)*P);
  for i=1:(nsym+1)*P
    st = 1 + (i-1)*Ts/P;
    en = st+Ts-1;
    f1_int(i) = sum(f1_dc(st:en));
    f2_int(i) = sum(f2_dc(st:en));
  end
  states.f1_int = f1_int;
  states.f2_int = f2_int;

  % fine timing estimation -----------------------------------------------

  % Non linearity has a spectral line at Rs, with a phase
  % related to the fine timing offset.  See:
  %   http://www.rowetel.com/blog/?p=3573 
  % We have sampled the integrator output at Fs=P samples/symbol, so
  % lets do a single point DFT at w = 2*pi*f/Fs = 2*pi*Rs/(P*Rs)

  Np = length(f1_int);
  w = 2*pi*(Rs)/(P*Rs);
  x = ((abs(f1_int)-abs(f2_int)).^2) * exp(-j*w*(0:Np-1))';
  norm_rx_timing = angle(x)/(2*pi);
  rx_timing = norm_rx_timing*P;

  states.x = x;
  states.rx_timing = rx_timing;
  states.norm_rx_timing = norm_rx_timing;

  % work out how many input samples we need on the next call. The aim
  % is to keep angle(x) away from the -pi/pi (+/- 0.5 fine timing
  % offset) discontinuity.  The side effect is to track sample clock
  % offsets

  next_nin = N;
  if norm_rx_timing > 0.25
     next_nin += Ts/2;
  end
  if norm_rx_timing < -0.25;
     next_nin -= Ts/2;
  end
  states.nin = next_nin;

  % Re-sample integrator outputs using fine timing estimate and linear interpolation

  low_sample = floor(rx_timing);
  fract = rx_timing - low_sample;
  high_sample = ceil(rx_timing);

  if verbose
    printf("rx_timing: %3.2f low_sample: %d high_sample: %d fract: %3.3f nin_next: %d\n", rx_timing, low_sample, high_sample, fract, next_nin);
  end

  f1_int_resample = zeros(1,nsym);
  f2_int_resample = zeros(1,nsym);
  rx_bits = zeros(1,nsym);
  for i=1:nsym
    st = i*P+1;
    f1_int_resample(i) = f1_int(st+low_sample)*(1-fract) + f1_int(st+high_sample)*fract;
    f2_int_resample(i) = f2_int(st+low_sample)*(1-fract) + f2_int(st+high_sample)*fract;
    %f1_int_resample(i) = f1_int(st+1);
    %f2_int_resample(i) = f2_int(st+1);
    rx_bits(i) = abs(f1_int_resample(i)) < abs(f2_int_resample(i));
  end

  states.f1_int_resample = f1_int_resample;
  states.f2_int_resample = f2_int_resample;

endfunction


% look for unique word and return index of first UW bit, or -1 if no UW found

function uw_start = find_uw(states, start_bit, rx_bits)
  uw = states.uw;

  mapped_rx_bits = 2*rx_bits - 1;
  found = 0;
  uw_start = -1;

  for i=start_bit:length(rx_bits) - length(uw)
    corr  = mapped_rx_bits(i:i+length(uw)-1) * uw';
    if !found && (corr >= states.uw_thresh)
      uw_start = i;
      found = 1;
    end
  end
endfunction


% Demo modem with simulated tx signal, add noise, channel impairments ----------------------

function run_sim
  frames = 100;
  EbNodB = 20;
  timing_offset = 0.0; % see resample() for clock offset below
  test_frame_mode = 4;
  fading = 1;          % modulates tx power at 5Hz with 20dB fade depth, 
                       % to simulate balloon rotating at end of mission
  df     = 1;          % tx tone freq drift in Hz/s

  more off
  rand('state',1); 
  randn('state',1);
  states = fsk_horus_init();
  N = states.N;
  P = states.P;
  Rs = states.Rs;
  nsym = states.nsym;
  Fs = states.Fs;
  states.df = df;

  EbNo = 10^(EbNodB/10);
  variance = states.Fs/(states.Rs*EbNo);

  % set up tx signal with payload bits based on test mode

  if test_frame_mode == 1
     % test frame of bits, which we repeat for convenience when BER testing
    test_frame = round(rand(1, states.nsym));
    tx_bits = [];
    for i=1:frames+1
      tx_bits = [tx_bits test_frame];
    end
  end
  if test_frame_mode == 2
    % random bits, just to make sure sync algs work on random data
    tx_bits = round(rand(1, states.nsym*(frames+1)));
  end
  if test_frame_mode == 3
    % ...10101... sequence
    tx_bits = zeros(1, states.nsym*(frames+1));
    tx_bits(1:2:length(tx_bits)) = 1;
  end

  if test_frame_mode == 4

    % load up a horus msg from disk and modulate that

    test_frame = load("horus_msg.txt");
    ltf = length(test_frame);
    ntest_frames = ceil((frames+1)*nsym/ltf);
    tx_bits = [];
    for i=1:ntest_frames
      tx_bits = [tx_bits test_frame];
    end
  end

  tx = fsk_horus_mod(states, tx_bits);

  tx = resample(tx, 1000, 1000); % simulated 1000ppm sample clock offset

  if fading
     ltx = length(tx);
     tx = tx .* (1.1 + cos(2*pi*5*(0:ltx-1)/Fs))'; % min amplitude 0.1, -20dB fade, max 3dB
  end

  noise = sqrt(variance/2)*(randn(length(tx),1) + j*randn(length(tx),1));
  rx    = tx + noise;

  % dump simulated rx file
  ftx=fopen("fsk_horus_rx.raw","wb"); rxg = rx*5000; fwrite(ftx, rxg, "short"); fclose(ftx);

  timing_offset_samples = round(timing_offset*states.Ts);
  st = 1 + timing_offset_samples;
  rx_bits_buf = zeros(1,2*nsym);
  Terrs = Tbits = 0;
  state = 0;
  x_log = [];
  norm_rx_timing_log = [];
  nerr_log = [];
  f1_int_resample_log = [];
  f2_int_resample_log = [];
  f1_log = f2_log = [];

  for f=1:frames

    % extract nin samples from input stream

    nin = states.nin;
    en = st + states.nin - 1;
    sf = rx(st:en);
    st += nin;

    % demodulate to stream of bits

    [rx_bits states] = fsk_horus_demod(states, sf);
    rx_bits_buf(1:nsym) = rx_bits_buf(nsym+1:2*nsym);
    rx_bits_buf(nsym+1:2*nsym) = rx_bits;
    norm_rx_timing_log = [norm_rx_timing_log states.norm_rx_timing];
    x_log = [x_log states.x];
    f1_int_resample_log = [f1_int_resample_log abs(states.f1_int_resample)];
    f2_int_resample_log = [f2_int_resample_log abs(states.f2_int_resample)];
    f1_log = [f1_log states.f1];
    f2_log = [f2_log states.f2];

    % frame sync based on min BER

    if test_frame_mode == 1
      nerrs_min = nsym;
      next_state = state;
      if state == 0
        for i=1:nsym
          error_positions = xor(rx_bits_buf(i:nsym+i-1), test_frame);
          nerrs = sum(error_positions);
          if nerrs < nerrs_min
            nerrs_min = nerrs;
            coarse_offset = i;
          end
        end
        if nerrs_min < 3
          next_state = 1;
          %printf("%d %d\n", coarse_offset, nerrs_min);
        end
      end

      if state == 1  
        error_positions = xor(rx_bits_buf(coarse_offset:coarse_offset+nsym-1), test_frame);
        nerrs = sum(error_positions);
        Terrs += nerrs;
        Tbits += nsym;
        nerr_log = [nerr_log nerrs];
      end

      state = next_state;
    end
  end

  if test_frame_mode == 1
    printf("frames: %d Tbits: %d Terrs: %d BER %3.2f\n", frames, Tbits, Terrs, Terrs/Tbits);
  end

  figure(1);
  plot(f1_int_resample_log,'+')
  hold on;
  plot(f2_int_resample_log,'g+')
  hold off;

  figure(2)
  clf
  m = max(abs(x_log));
  plot(x_log,'+')
  axis([-m m -m m])
  title('fine timing metric')

  figure(3)
  clf
  subplot(211)
  plot(norm_rx_timing_log);
  axis([1 frames -1 1])
  title('norm fine timing')
  subplot(212)
  plot(nerr_log)
  title('num bit errors each frame')

  figure(4)
  clf
  plot(real(rx(1:Fs)))
  title('rx signal at demod input')

  figure(5)
  clf
  plot(f1_log)
  hold on;
  plot(f2_log,'g');
  hold off;
  title('tone frequencies')
  axis([1 frames 0 Fs/2])
endfunction


% demodulate a file of 8kHz 16bit short samples --------------------------------

function rx_bits_log = demod_file(filename)
  rx = load_raw(filename);
  more off;
  states = fsk_horus_init();
  N = states.N;
  P = states.P;
  Rs = states.Rs;
  nsym = states.nsym;

  frames = floor(length(rx)/N);
  st = 1;
  rx_bits_log = [];
  norm_rx_timing_log = [];
  f1_int_resample_log = [];
  f2_int_resample_log = [];

  % First extract raw bits from samples ------------------------------------------------------

  printf("demod of raw bits....\n");

  for f=1:frames

    % extract nin samples from input stream

    nin = states.nin;
    en = st + states.nin - 1;
    sf = rx(st:en);
    st += nin;

    % demodulate to stream of bits

    [rx_bits states] = fsk_horus_demod(states, sf);
    rx_bits_log = [rx_bits_log rx_bits];
    norm_rx_timing_log = [norm_rx_timing_log states.norm_rx_timing];
    f1_int_resample_log = [f1_int_resample_log abs(states.f1_int_resample)];
    f2_int_resample_log = [f2_int_resample_log abs(states.f2_int_resample)];
  end

  printf("plotting...\n");

  figure(1);
  plot(f1_int_resample_log,'+')
  hold on;
  plot(f2_int_resample_log,'g+')
  hold off;

  figure(2)
  clf
  plot(norm_rx_timing_log)
  axis([1 frames -0.5 0.5])
  title('norm fine timing')
  grid
  
  printf("frame sync and data extraction...\n");

  % Now perform frame sync and extract ASCII text -------------------------------------------

  % use UWs to delimit start and end of data packets

  bit = 1;
  nbits = length(rx_bits_log);
  uw_loc = find_uw(states, bit, rx_bits_log);
  nfield = states.nfield;
  npad = states.npad;

  while (uw_loc != -1)

    st = uw_loc;
    bit = uw_loc + length(states.uw);
    uw_loc = find_uw(states, bit, rx_bits_log);

    if uw_loc != -1
      % Now start picking out 7 bit ascii chars from frame.  It has some
      % structure so we can guess where fields are.  I hope We don't get
      % RS232 idle bits stuck into it anywhere, ie "bit fields" don't
      % change dynamically.

      % dump msg bits so we can use them as a test signal
      %msg = rx_bits_log(st:uw_loc-1);
      %save -ascii horus_msg.txt msg

      str = [];
      st += length(states.uw);  % first bit of first char
      for i=st:nfield+npad:uw_loc
        field = rx_bits_log(i:i+nfield-1);
        ch_dec = field * (2.^(0:nfield-1))';
        % filter out unlikely characters that bit errors may introduce, and ignore \n
        if (ch_dec > 31) && (ch_dec < 91)
          str = [str char(ch_dec)];
        else 
          str = [str char(32)]; % space is "non sure"
        end
        %printf("i: %d ch_dec: %d ch: %c\n", i, ch_dec, ch_dec);
      end
      printf("%s\n", str);
    end
   
  endwhile
 
endfunction


% run test functions from here during development

%run_sim
%rx_bits = demod_file("~/Desktop/vk5arg-3.wav");
%rx_bits = demod_file("~/Desktop/fsk_horus_10dB_1000ppm.wav");
%rx_bits = demod_file("~/Desktop/fsk_horus_6dB_0ppm.wav");
%rx_bits = demod_file("fsk_horus_rx.raw");
%rx_bits = demod_file("~/Desktop/fsk_horus_20dB_0ppm_20dBfade.wav");

% streaming version: take a buffer of say 1 sec.  demo to bits.  Add to buffer
% of bits.  When two UWs found, demod and dump txt.  Shift buffer back that far.
% how long to make buffer?  How to dunmp diagnostics?  printf f1, f2, clock offset est
