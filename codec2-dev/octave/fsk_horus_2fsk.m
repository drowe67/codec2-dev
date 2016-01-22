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
% [X] test over range of f1/f2, shifts, timing offsets, clock offsets, Eb/No
%     [X] +/- 1000ppm clock offset OK at Eb/No = 10dB, starts to lose it at 8dB
%     [X] tone freq est starts to lose it at 8dB in awgn.  Maybe wider window?
% [ ] low snr detection of $$$$$$
%     + we might be able to pick up a "ping" at very low SNRs to help find baloon on ground
% [ ] streaming, indicator of audio freq, i.e. speaker output?

% horus binary:
% [ ] BER estimate/found/corrected

1;

function states = fsk_horus_init(Fs,Rs)
  states.Ndft = 2.^ceil(log2(Fs)); % find nearest power of 2 for effcient FFT
  states.Fs = Fs;
  N = states.N = Fs;             % processing buffer size, nice big window for f1,f2 estimation
  states.Rs = Rs;
  Ts = states.Ts = Fs/Rs;
  assert(Ts == floor(Ts), "Fs/Rs must be an integer");
  states.nsym = N/Ts;            % number of symbols in one procesing frame
  Nmem = states.Nmem  = N+2*Ts;  % two symbol memory in down converted signals to allow for timing adj
  states.f1_dc = zeros(1,Nmem);
  states.f2_dc = zeros(1,Nmem);
  states.P = 8;                  % oversample rate out of filter
  assert(Ts/states.P == floor(Ts/states.P), "Ts/P must be an integer");

  states.nin = N;                % can be N +/- Ts/P samples to adjust for sample clock offsets
  states.verbose = 0;
  states.phi1 = 0;               % keep down converter osc phase continuous
  states.phi2 = 0;

  printf("Fs: %d Rs: %d Ts: %d nsym: %d\n", states.Fs, states.Rs, states.Ts, states.nsym);

  % BER stats 

  states.ber_state = 0;
  states.Tbits = 0;
  states.Terrs = 0;
  states.nerr_log = 0;

  states.df = 0;
  states.f1 = 0;
  states.f2 = 0;
  states.norm_rx_timing = 0;
  states.ppm = 0;
  states.prev_pkt = [];

  % protocol specific states

  states.rtty = fsk_horus_init_rtty_uw(states);
  states.binary = fsk_horus_init_binary_uw;
endfunction


% init rtty protocol specifc states

function rtty = fsk_horus_init_rtty_uw(states)
  % Generate unque word that correlates against the ASCII "$$$$$" that
  % is at the start of each frame.

  dollar_bits = fliplr([0 1 0 0 1 0 0]);
  mapped_db = 2*dollar_bits - 1;
  sync_bits = [1 1 0];
  mapped_sb = 2*sync_bits - 1;
  %mapped_sb = [ 0 0 0 ];

  mapped = [mapped_db mapped_sb];
  npad = rtty.npad = 3;     % one start and two stop bits between 7 bit ascii chars
  nfield = rtty.nfield = 7; % length of ascii character field

  rtty.uw = [mapped mapped mapped mapped mapped];

  rtty.uw_thresh = length(rtty.uw) - 8; % allow a few bit errors when looking for UW

  rtty.max_packet_len = 1000;
endfunction



function binary = fsk_horus_init_binary_uw
  % Generate 16 bit "$$" unique word that is at the front of every horus binary
  % packet

  dollar_bits = [0 0 1 0 0 1 0 0];
  mapped_db = 2*dollar_bits - 1;

  binary.uw = [mapped_db mapped_db];
  binary.uw_thresh = length(binary.uw);   % no bit errors when looking for UW

  binary.max_packet_len = 400;
endfunction


% test modulator function

function tx  = fsk_horus_mod(states, tx_bits)
    tx = zeros(states.Ts*length(tx_bits),1);
    tx_phase = 0;
    Ts = states.Ts;
    Fs = states.Fs;
    f1 = states.f1_tx;  f2 = states.f2_tx;
    df = states.df; % tone freq change in Hz/s
    dA = states.dA;

    for i=1:length(tx_bits)
      if tx_bits(i) == 0
        tx_phase_vec = tx_phase + (1:Ts)*2*pi*f1/Fs;
        tx((i-1)*Ts+1:i*Ts) = dA*2.0*cos(tx_phase_vec);
      else
        tx_phase_vec = tx_phase + (1:Ts)*2*pi*f2/Fs;
        tx((i-1)*Ts+1:i*Ts) = 2.0*cos(tx_phase_vec);
      end
      tx_phase = tx_phase_vec(Ts) - floor(tx_phase_vec(Ts)/(2*pi))*2*pi;
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

  % zero out region 100Hz either side of max so we can find second highest peak

  Sf2 = Sf;
  st = m1_index - floor(100*Ndft/Fs);
  if st < 1
    st = 1;
  end
  en = m1_index + floor(100*Ndft/Fs);
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

  if bitand(verbose,0x1)
    printf("centre: %4.0f shift: %4.0f twist: %3.1f dB\n", (f2+f1)/2, f2-f1, twist);
  end
  if bitand(verbose,0x8)
    printf("f1: %4.0f Hz f2: %4.0f Hz a1: %f a2: %f\n", f1, f2, 2.0*abs(m1)/Ndft, 2.0*abs(m2)/Ndft);
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
  prev_norm_rx_timing = states.norm_rx_timing;
  states.norm_rx_timing = norm_rx_timing;

  % estimate sample clock offset in ppm
  % d_norm_timing is fraction of symbol period shift over nsym symbols

  d_norm_rx_timing = norm_rx_timing - prev_norm_rx_timing;

  % filter out big jumps due to nin changes

  if abs(d_norm_rx_timing) < 0.2
    appm = 1E6*d_norm_rx_timing/nsym;
    states.ppm = 0.9*states.ppm + 0.1*appm;
  end

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

  if bitand(verbose,0x2)
    printf("rx_timing: %3.2f low_sample: %d high_sample: %d fract: %3.3f nin_next: %d\n", rx_timing, low_sample, high_sample, fract, next_nin);
  end

  f1_int_resample = zeros(1,nsym);
  f2_int_resample = zeros(1,nsym);
  rx_bits = zeros(1,nsym);
  rx_bits_sd = zeros(1,nsym);
  for i=1:nsym
    st = i*P+1;
    f1_int_resample(i) = f1_int(st+low_sample)*(1-fract) + f1_int(st+high_sample)*fract;
    f2_int_resample(i) = f2_int(st+low_sample)*(1-fract) + f2_int(st+high_sample)*fract;
    %f1_int_resample(i) = f1_int(st+1);
    %f2_int_resample(i) = f2_int(st+1);
    rx_bits(i) = abs(f2_int_resample(i)) > abs(f1_int_resample(i));
    rx_bits_sd(i) = abs(f2_int_resample(i)) - abs(f1_int_resample(i));
 end

  states.f1_int_resample = f1_int_resample;
  states.f2_int_resample = f2_int_resample;
  states.rx_bits_sd = rx_bits_sd;

  % Eb/No estimation

  x = abs(abs(f1_int_resample) - abs(f2_int_resample));
  states.EbNodB = 20*log10(1E-6+mean(x)/(1E-6+std(x)));
endfunction


% Look for unique word and return index of first UW bit, or -1 if no
% UW found Sometimes there may be several matches, returns the
% position of the best match to UW.

function uw_start = find_uw(states, start_bit, rx_bits)
  uw = states.uw;

  mapped_rx_bits = 2*rx_bits - 1;
  best_corr = 0;
  uw_start = -1;

  for i=start_bit:length(rx_bits) - length(uw)
    corr  = mapped_rx_bits(i:i+length(uw)-1) * uw';
    if (corr >= states.uw_thresh) && (corr > best_corr)
      uw_start = i;
      best_corr = corr;
    end
  end
endfunction


% Extract ASCII string from a Horus frame of bits

function [str crc_ok] = extract_ascii(states, rx_bits_buf, uw_loc)
  nfield = states.nfield;
  npad = states.npad;

  str = []; str_dec = []; nstr = 0; ptx_crc = 1; rx_crc = "";
  endpacket = 0;

  st = uw_loc + length(states.uw);  % first bit of first char
  en = uw_loc + states.max_packet_len - nfield;
  %printf("\nst: %d en: %d len: %d\n", st, en, length(rx_bits_buf));

  for i=st:nfield+npad:en
    field = rx_bits_buf(i:i+nfield-1);
    ch_dec = field * (2.^(0:nfield-1))';

    % filter out unlikely characters that bit errors may introduce, and ignore \n

    if (ch_dec > 31) && (ch_dec < 91)
      str = [str char(ch_dec)];
    else 
      str = [str char(32)]; % space is "not sure"
    end
    nstr++;

    % build up array for CRC16 check

    if !endpacket && (ch_dec == 42)
      endpacket = 1; 
      rx_crc = crc16(str_dec);      % found a '*' so that's the end of the string for CRC calculations
      ptx_crc = nstr+1;             % this is where the transmit CRC starts
    end
    if !endpacket
      str_dec = [str_dec ch_dec];
    end
  end

  if (ptx_crc+3) <= length(str)
    tx_crc = str(ptx_crc:ptx_crc+3);
    crc_ok = strcmp(tx_crc, rx_crc);
  else
    crc_ok = 0;
  end

  str = str(1:ptx_crc-2);

endfunction


% Use soft decision information to find bits most likely in error.  I think
% this is some form of maximum likelihood decoding.

function [str crc_ok rx_bits_log_flipped] = sd_bit_flipping(states, rx_bits_log, rx_bits_sd_log, st, en);

  % force algorithm to ignore rs232 sync bits by marking them as "very likely", they have
  % no input to crc algorithm

  nfield = states.nfield;
  npad = states.npad;
  for i=st:nfield+npad:en
    rx_bits_sd_log(i+nfield:i+nfield+npad-1) = 1E6;
  end

  % make a list of bits with smallest soft decn values

  [dodgy_bits_mag dodgy_bits_index] = sort(abs(rx_bits_sd_log(st+length(states.uw):en)));
  dodgy_bits_index += length(states.uw) + st - 1;
  nbits = 6;
  ntries = 2^nbits;
  str = "";
  crc_ok = 0;
  
  % try various combinations of these bits

  for i=1:ntries-1
    error_mask = zeros(1, length(rx_bits_log));
    for b=1:nbits
      x = bitget(i,b);
      bit_to_flip = dodgy_bits_index(b);
      error_mask(bit_to_flip) = x;
      %printf("st: %d i: %d b: %d x: %d index: %d\n", st, i,b,x,bit_to_flip);
    end
    rx_bits_log_flipped = xor(rx_bits_log, error_mask);
    [str_flipped crc_ok_flipped] = extract_ascii(states, rx_bits_log_flipped, st);
    if crc_ok_flipped
      %printf("Yayy we fixed a packet by flipping with pattern %d\n", i);
      str = str_flipped;
      crc_ok = crc_ok_flipped;
    end
  end
endfunction


% Extract as many ASCII packets as we can from a great big buffer of bits

function extract_and_print_rtty_packets(states, rx_bits_log, rx_bits_sd_log)

  % use UWs to delimit start and end of data packets

  bit = 1;
  nbits = length(rx_bits_log);
  nfield = states.rtty.nfield;
  npad = states.rtty.npad;

  uw_loc = find_uw(states.rtty, bit, rx_bits_log);

  while (uw_loc != -1)

    if (uw_loc + states.rtty.max_packet_len) < nbits
      % Now start picking out 7 bit ascii chars from frame.  It has some
      % structure so we can guess where fields are.  I hope we don't get
      % RS232 idle bits stuck into it anywhere, ie "bit fields" don't
      % change dynamically.

      % dump msg bits so we can use them as a test signal
      %msg = rx_bits_log(st:uw_loc-1);
      %save -ascii horus_msg.txt msg

      % simulate bit error for testing
      %rx_bits_log(st+200) = xor(rx_bits_log(st+100),1);
      %rx_bits_sd_log(st+100) = 0;
      
      [str crc_ok] = extract_ascii(states.rtty, rx_bits_log, uw_loc);

      if crc_ok == 0
        [str_flipped crc_flipped_ok rx_bits_log] = sd_bit_flipping(states.rtty, rx_bits_log, rx_bits_sd_log, uw_loc, uw_loc+states.rtty.max_packet_len); 
        if crc_flipped_ok
          str = sprintf("%s fixed", str_flipped);
        end
      end

      % update memory of previous packet, we use this to guess where errors may be
      if crc_ok || crc_flipped_ok
        states.prev_pkt = rx_bits_log(uw_loc+length(states.rtty.uw):uw_loc+states.rtty.max_packet_len);
      end
      if crc_ok
        str = sprintf("%s CRC OK", str);
      else
        str = sprintf("%s CRC BAD", str);
      end
      printf("%s\n", str);
    end

    % look for next packet

    bit = uw_loc + length(states.rtty.uw);
    uw_loc = find_uw(states.rtty, bit, rx_bits_log);

  endwhile
endfunction
 

% Extract as many ASCII packets as we can from a great big buffer of bits,
% and send them to the C decoder for FEC decoding.
% horus_l2 can be compiled a bunch of different ways.  You need to
% compile with:
%   codec2-dev/src$ gcc horus_l2.c -o horus_l2 -Wall -DDEC_RX_BITS -DHORUS_L2_RX

function extract_and_decode_binary_packets(states, rx_bits_log)

  % use UWs to delimit start and end of data packets

  bit = 1;
  nbits = length(rx_bits_log);

  uw_loc = find_uw(states.binary, bit, rx_bits_log);

  while (uw_loc != -1)

    if (uw_loc+states.binary.max_packet_len) < nbits
      %printf("st: %d uw_loc: %d\n", st, uw_loc);

      % OK we have a packet delimited by two UWs.  Lets convert the bit
      % stream into bytes and save for decoding

      pin = uw_loc;    
      for i=1:45
        rx_bytes(i) = rx_bits_log(pin:pin+7) * (2.^(7:-1:0))';
        pin += 8;
        %printf("%d 0x%02x\n", i, rx_bytes(i));
      end

      f=fopen("horus_rx_bits_binary.txt","wt");
      fwrite(f, rx_bytes, "uchar");
      fclose(f);

      system("../src/horus_l2");  % compile instructions above
    end

    bit = uw_loc + length(states.binary.uw);
    uw_loc = find_uw(states.binary, bit, rx_bits_log);
   
  endwhile
endfunction
 

% BER counter and test frame sync logic

function states = ber_counter(states, test_frame, rx_bits_buf)
  nsym = states.nsym;
  state = states.ber_state;
  next_state = state;

  if state == 0

    % try to sync up with test frame

    nerrs_min = nsym;
    for i=1:nsym
      error_positions = xor(rx_bits_buf(i:nsym+i-1), test_frame);
      nerrs = sum(error_positions);
      if nerrs < nerrs_min
        nerrs_min = nerrs;
        states.coarse_offset = i;
      end
    end
    if nerrs_min/nsym < 0.05 
      next_state = 1;
    end
    if bitand(states.verbose,0x4)
      printf("coarse offset: %d nerrs_min: %d next_state: %d\n", states.coarse_offset, nerrs_min, next_state);
    end
  end

  if state == 1  

    % we're synced up, lets measure bit errors

    error_positions = xor(rx_bits_buf(states.coarse_offset:states.coarse_offset+nsym-1), test_frame);
    nerrs = sum(error_positions);
    if nerrs/nsym > 0.1
      next_state = 0;
    else
      states.Terrs += nerrs;
      states.Tbits += nsym;
      states.nerr_log = [states.nerr_log nerrs];
    end
  end

  states.ber_state = next_state;
endfunction


% simulation of tx and rx side, add noise, channel impairments ----------------------

function run_sim(test_frame_mode)
  frames = 60;
  EbNodB = 10;
  timing_offset = 0.0; % see resample() for clock offset below
  fading = 0;          % modulates tx power at 2Hz with 20dB fade depth, 
                       % to simulate balloon rotating at end of mission
  df     = 0;          % tx tone freq drift in Hz/s
  dA     = 1;          % amplitude imbalance of tones (note this affects Eb so not a gd idea)

  more off
  rand('state',1); 
  randn('state',1);

  % ----------------------------------------------------------------------

  % sm2000 config ------------------------
  %states = fsk_horus_init(96000, 1200);
  %states.f1_tx = 4000;
  %states.f2_tx = 5200;

  if test_frame_mode == 4
    % horus rtty config ---------------------
    states = fsk_horus_init(8000, 100);
    states.f1_tx = 1200;
    states.f2_tx = 1600;
    states.tx_bits_file = "horus_tx_bits_rtty.txt"; % Octave file of bits we FSK modulate
  end
                               
  if test_frame_mode == 5
    % horus binary config ---------------------
    states = fsk_horus_init(8000, 100);
    states.f1_tx = 1200;
    states.f2_tx = 1600;
    states.tx_bits_file = "horus_tx_bits_binary.txt"; % Octave file of bits we FSK modulate
  end

  % ----------------------------------------------------------------------

  states.verbose = 0x1;
  N = states.N;
  P = states.P;
  Rs = states.Rs;
  nsym = states.nsym;
  Fs = states.Fs;
  states.df = df;
  states.dA = dA;

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
 
  if (test_frame_mode == 4) || (test_frame_mode == 5)

    % load up a horus msg from disk and modulate that

    test_frame = load(states.tx_bits_file);
    ltf = length(test_frame);
    ntest_frames = ceil((frames+1)*nsym/ltf);
    tx_bits = [];
    for i=1:ntest_frames
      tx_bits = [tx_bits test_frame];
    end
  end

  tx = fsk_horus_mod(states, tx_bits);

  %tx = resample(tx, 1000, 1001); % simulated 1000ppm sample clock offset

  if fading
     ltx = length(tx);
     tx = tx .* (1.1 + cos(2*pi*2*(0:ltx-1)/Fs))'; % min amplitude 0.1, -20dB fade, max 3dB
  end

  noise = sqrt(variance)*randn(length(tx),1);
  rx    = tx + noise;
  %rx = real(rx);
  %b1 = fir2(100, [0 4000 5200 48000]/48000, [1 1 0.5 0.5]);
  %rx = filter(b1,1,rx);
  %[b a] = cheby2(6,40,[3000 6000]/(Fs/2));
  %rx = filter(b,a,rx);
  %rx = sign(rx);
  %rx(find (rx > 1)) = 1;
  %rx(find (rx < -1)) = -1;

  % dump simulated rx file
  ftx=fopen("fsk_horus_100bd_binary.raw","wb"); rxg = rx*1000; fwrite(ftx, rxg, "short"); fclose(ftx);

  timing_offset_samples = round(timing_offset*states.Ts);
  st = 1 + timing_offset_samples;
  rx_bits_buf = zeros(1,2*nsym);
  x_log = [];
  norm_rx_timing_log = [];
  f1_int_resample_log = [];
  f2_int_resample_log = [];
  f1_log = f2_log = [];
  EbNodB_log = [];
  rx_bits_log = [];
  rx_bits_sd_log = [];

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
    rx_bits_log = [rx_bits_log rx_bits];
    rx_bits_sd_log = [rx_bits_sd_log states.rx_bits_sd];

    norm_rx_timing_log = [norm_rx_timing_log states.norm_rx_timing];
    x_log = [x_log states.x];
    f1_int_resample_log = [f1_int_resample_log abs(states.f1_int_resample)];
    f2_int_resample_log = [f2_int_resample_log abs(states.f2_int_resample)];
    f1_log = [f1_log states.f1];
    f2_log = [f2_log states.f2];
    EbNodB_log = [EbNodB_log states.EbNodB];

    if test_frame_mode == 1
       states = ber_counter(states, test_frame, rx_bits_buf);
    end
  end

  if test_frame_mode == 1
    printf("frames: %d Tbits: %d Terrs: %d BER %4.3f\n", frames, states.Tbits,states. Terrs, states.Terrs/states.Tbits);
  end

  if test_frame_mode == 4
    extract_and_print_rtty_packets(states, rx_bits_log, rx_bits_sd_log)
  end

  if test_frame_mode == 5
    extract_and_decode_binary_packets(states, rx_bits_log);
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
  plot(states.nerr_log)
  title('num bit errors each frame')

  figure(4)
  clf
  subplot(211)
  plot(real(rx(1:Fs)))
  title('rx signal at demod input')
  subplot(212)
  plot(abs(fft(rx(1:Fs))))

  figure(5)
  clf
  plot(f1_log)
  hold on;
  plot(f2_log,'g');
  hold off;
  title('tone frequencies')
  axis([1 frames 0 Fs/2])

  figure(6)
  clf
  plot(EbNodB_log);
  title('Eb/No estimate')

endfunction


% demodulate a file of 8kHz 16bit short samples --------------------------------

function rx_bits_log = demod_file(filename, test_frame_mode, noplot)
  fin = fopen(filename,"rb"); 
  more off;

  %states = fsk_horus_init(96000, 1200);

  if test_frame_mode == 4
    % horus rtty config ---------------------
    states = fsk_horus_init(8000, 100);
    uwstates = fsk_horus_init_rtty_uw(states);
  end
                               
  if test_frame_mode == 5
    % horus binary config ---------------------
    states = fsk_horus_init(8000, 100);
    uwstates = fsk_horus_init_binary_uw;
  end

  states.verbose = 0x1 + 0x8;

  N = states.N;
  P = states.P;
  Rs = states.Rs;
  nsym = states.nsym;
  rand('state',1); 
  test_frame = round(rand(1, states.nsym));

  frames = 0;
  rx = [];
  rx_bits_log = [];
  rx_bits_sd_log = [];
  norm_rx_timing_log = [];
  f1_int_resample_log = [];
  f2_int_resample_log = [];
  EbNodB_log = [];
  ppm_log = [];
  rx_bits_buf = zeros(1,2*nsym);

  % First extract raw bits from samples ------------------------------------------------------

  printf("demod of raw bits....\n");

  finished = 0;
  while (finished == 0)

    % hit any key to finish (useful for real time streaming)

    %x = kbhit(1);
    %if length(x)
    %  finished = 1;
    %end

    % extract nin samples from input stream

    nin = states.nin;
    [sf count] = fread(fin, nin, "short");
    rx = [rx; sf];

    if count == nin
      frames++;

      % demodulate to stream of bits

      [rx_bits states] = fsk_horus_demod(states, sf);
      rx_bits_buf(1:nsym) = rx_bits_buf(nsym+1:2*nsym);
      rx_bits_buf(nsym+1:2*nsym) = rx_bits; % xor(rx_bits,ones(1,nsym));
      rx_bits_log = [rx_bits_log rx_bits];
      rx_bits_sd_log = [rx_bits_sd_log states.rx_bits_sd];
      norm_rx_timing_log = [norm_rx_timing_log states.norm_rx_timing];
      f1_int_resample_log = [f1_int_resample_log abs(states.f1_int_resample)];
      f2_int_resample_log = [f2_int_resample_log abs(states.f2_int_resample)];
      EbNodB_log = [EbNodB_log states.EbNodB];
      ppm_log = [ppm_log states.ppm];

      if test_frame_mode == 1
        states = ber_counter(states, test_frame, rx_bits_buf);
        if states.ber_state == 1
          states.verbose = 0;
        end
      end
    else
      finished = 1;
    end
  end
  fclose(fin);

  if exist("noplot") == 0
    printf("plotting...\n");

    figure(1);
    plot(f1_int_resample_log,'+')
    hold on;
    plot(f2_int_resample_log,'g+')
    hold off;

    figure(2)
    clf
    subplot(211)
    plot(norm_rx_timing_log)
    axis([1 frames -0.5 0.5])
    title('norm fine timing')
    grid
    subplot(212)
    plot(states.nerr_log)
    title('num bit errors each frame')
 
    figure(3)
    clf
    plot(EbNodB_log);
    title('Eb/No estimate')

    figure(4)
    clf
    subplot(211)
    plot(rx(1:states.Fs));
    title('input signal to demod (1 sec)')
    xlabel('Time (samples)');
    axis([1 states.Fs -35000 35000])

    % normalise spectrum to 0dB full scale witha 32767 sine wave input

    subplot(212)
    RxdBFS = 20*log10(abs(fft(rx(1:states.Fs)))) - 20*log10((states.Fs/2)*32767);
    plot(RxdBFS)
    axis([1 states.Fs/2 -80 0])
    xlabel('Frequency (Hz)');

    figure(5);
    clf
    plot(ppm_log)
    title('Sample clock (baud rate) offset in PPM');
  end

  if test_frame_mode == 1
    printf("frames: %d Tbits: %d Terrs: %d BER %4.3f EbNo: %3.2f\n", frames, states.Tbits,states. Terrs, states.Terrs/states.Tbits, mean(EbNodB_log));
  end

  % we can decode both protocols at the same time

  if (test_frame_mode == 4) || (test_frame_mode == 5)
    extract_and_print_rtty_packets(states, rx_bits_log, rx_bits_sd_log)
    extract_and_decode_binary_packets(states, rx_bits_log);
  end
endfunction


% run test functions from here during development

if exist("fsk_horus_as_a_lib") == 0
  %run_sim(5);
  %rx_bits = demod_file("horus.raw",4);
  %rx_bits = demod_file("fsk_horus_100bd_binary.raw",5);
  rx_bits = demod_file("~/Desktop/phorus_binary_ascii.wav",4);
  %rx_bits = demod_file("~/Desktop/horus_rtty_binary.wav",4);
  %rx_bits = demod_file("t.raw",5);
  %rx_bits = demod_file("~/Desktop/fsk_horus_10dB_1000ppm.wav",4);
  %rx_bits = demod_file("~/Desktop/fsk_horus_6dB_0ppm.wav",4);
  %rx_bits = demod_file("test.raw",1,1);
  %rx_bits = demod_file("/dev/ttyACM0",1);
  %rx_bits = demod_file("fsk_horus_rx_1200_96k.raw",1);
  %rx_bits = demod_file("mp.raw",4);
  %rx_bits = demod_file("~/Desktop/launchbox_v2_landing_8KHz_final.wav",4);
end
