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

function states = fsk_horus_init(Fs,Rs,M=2)
  assert((M==2) || (M==4), "Only M=2 and M=4 FSK supported");
  states.M = M;                    
  states.bitspersymbol = log2(M);
  states.Fs = Fs;
  N = states.N = Fs;                % processing buffer size, nice big window for timing est
  %states.Ndft = 2.^ceil(log2(N));  % find nearest power of 2 for efficient FFT
  states.Ndft = 1024;               % find nearest power of 2 for efficient FFT
  states.Rs = Rs;
  Ts = states.Ts = Fs/Rs;
  assert(Ts == floor(Ts), "Fs/Rs must be an integer");
  states.nsym = N/Ts;            % number of symbols in one processing frame
  states.nbit = states.nsym*states.bitspersymbol; % number of bits per processing frame

  Nmem = states.Nmem  = N+2*Ts;  % two symbol memory in down converted signals to allow for timing adj

  states.Sf = zeros(states.Ndft/2,1); % currentmemory of dft mag samples
  states.f_dc = zeros(M,Nmem);
  states.P = 8;                  % oversample rate out of filter
  assert(Ts/states.P == floor(Ts/states.P), "Ts/P must be an integer");

  states.nin = N;                % can be N +/- Ts/P samples to adjust for sample clock offsets
  states.verbose = 0;
  states.phi = zeros(1, M);       % keep down converter osc phase continuous

  %printf("M: %d Fs: %d Rs: %d Ts: %d nsym: %d nbit: %d\n", states.M, states.Fs, states.Rs, states.Ts, states.nsym, states.nbit);

  % BER stats 

  states.ber_state = 0;
  states.Tbits = 0;
  states.Terrs = 0;
  states.nerr_log = 0;

  % extra simulation parameters

  states.tx_real = 1;
  states.dA(1:M) = 1;
  states.df(1:M) = 0;
  states.f(1:M) = 0;
  states.norm_rx_timing = 0;
  states.ppm = 0;
  states.prev_pkt = [];
 
  % protocol specific states

  states.rtty = fsk_horus_init_rtty_uw(states);
  states.binary = fsk_horus_init_binary_uw;

  % Freq. estimator limits - keep these narrow to stop errors with low SNR 4FSK

  states.fest_fmin = 800;
  states.fest_fmax = 2500;
  states.fest_min_spacing = 200;
endfunction


% Alternative init function, useful for high speed (non telemetry) modems
%   Allows fine grained control of decimation P
%   Small, processing window nsym rather than nsym=Fs (1 second window)
%   Wider freq est limits

function states = fsk_horus_init_hbr(Fs,P,Rs,M=2,nsym=48)
  assert((M==2) || (M==4), "Only M=2 and M=4 FSK supported");
    
  states.M = M;                    
  states.bitspersymbol = log2(M);
  states.Fs = Fs;
  states.Rs = Rs;
  Ts = states.Ts = Fs/Rs;
  assert(Ts == floor(Ts), "Fs/Rs must be an integer");
  N = states.N = Ts*nsym;        % processing buffer nsym wide
  states.nsym = N/Ts;            % number of symbols in one processing frame
  states.nbit = states.nsym*states.bitspersymbol; % number of bits per processing frame

  states.Ndft = (2.^ceil(log2(N)))/2;  % find nearest power of 2 for efficient FFT

  Nmem = states.Nmem  = N+2*Ts;  % two symbol memory in down converted signals to allow for timing adj

  states.Sf = zeros(states.Ndft/2,1); % currentmemory of dft mag samples
  states.f_dc = zeros(M,Nmem);
  states.P = P;                  % oversample rate out of filter
  assert(Ts/states.P == floor(Ts/states.P), "Ts/P must be an integer");

  states.nin = N;                % can be N +/- Ts/P samples to adjust for sample clock offsets
  states.verbose = 0;
  states.phi = zeros(1, M);      % keep down converter osc phase continuous

  %printf("M: %d Fs: %d Rs: %d Ts: %d nsym: %d nbit: %d\n", states.M, states.Fs, states.Rs, states.Ts, states.nsym, states.nbit);

  % Freq estimator limits

  states.fest_fmax = (Fs/2)-Rs;
  states.fest_fmin = Rs/2;
  states.fest_min_spacing = 2*(Rs-(Rs/5));

  % BER stats 

  states.ber_state = 0;
  states.Tbits = 0;
  states.Terrs = 0;
  states.nerr_log = 0;

  states.tx_real = 1;
  states.dA(1:M) = 1;
  states.df(1:M) = 0;
  states.f(1:M) = 0;
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
  % $ -> 36 decimal -> 0 1 0 0 1 0 0 binary 

  dollar_bits = fliplr([0 1 0 0 1 0 0]);
  mapped_db = 2*dollar_bits - 1;
  sync_bits = [1 1 0];
  mapped_sb = 2*sync_bits - 1;
  %mapped_sb = [ 0 0 0 ];

  mapped = [mapped_db mapped_sb];
  npad = rtty.npad = 3;     % one start and two stop bits between 7 bit ascii chars
  nfield = rtty.nfield = 7; % length of ascii character field

  rtty.uw = [mapped mapped mapped mapped mapped];
  rtty.uw_thresh = length(rtty.uw) - 2; % allow a few bit errors when looking for UW
  rtty.max_packet_len = 1000;
endfunction



function binary = fsk_horus_init_binary_uw
  % Generate 16 bit "$$" unique word that is at the front of every horus binary
  % packet

  dollar_bits = [0 0 1 0 0 1 0 0];
  mapped_db = 2*dollar_bits - 1;

  binary.uw = [mapped_db mapped_db];
  binary.uw_thresh = length(binary.uw)-2;   % no bit errors when looking for UW

  binary.max_packet_len = 360;
endfunction


% test modulator function

function tx  = fsk_horus_mod(states, tx_bits)

    M  = states.M;
    Ts = states.Ts;
    Fs = states.Fs;
    ftx  = states.ftx;
    df = states.df; % tone freq change in Hz/s
    dA = states.dA; % amplitude of each tone

    num_bits = length(tx_bits);
    num_symbols = num_bits/states.bitspersymbol;
    tx = zeros(states.Ts*num_symbols,1);
    tx_phase = 0;
    s = 1;

    for i=1:states.bitspersymbol:num_bits

      % map bits to tone number

      if M == 2
        tone = tx_bits(i) + 1;
      else
        tone = (tx_bits(i:i+1) * [2 1]') + 1;
      end
 
      tx_phase_vec = tx_phase + (1:Ts)*2*pi*ftx(tone)/Fs;
      tx_phase = tx_phase_vec(Ts) - floor(tx_phase_vec(Ts)/(2*pi))*2*pi;
      if states.tx_real
        tx((s-1)*Ts+1:s*Ts) = dA(tone)*2.0*cos(tx_phase_vec);
      else
        tx((s-1)*Ts+1:s*Ts) = dA(tone)*exp(j*tx_phase_vec);
      end
      s++;

      % freq drift

      ftx += df*Ts/Fs;
    end
    states.ftx = ftx;
endfunction


% Estimate the frequency of the FSK tones.  In some applications (such
% as balloon telemtry) these may not be well controlled by the
% transmitter, so we have to try to estimate them.

function states = est_freq(states, sf, ntones)
  N = states.N;
  Ndft = states.Ndft;
  Fs = states.Fs;
  
  % This assumption is OK for balloon telemetry but may not be true in
  % general

  min_tone_spacing = states.fest_min_spacing;
  
  % set some limits to search range, which will mean some manual re-tuning

  fmin = states.fest_fmin;
  fmax = states.fest_fmax;
  st = floor(fmin*Ndft/Fs);
  en = floor(fmax*Ndft/Fs);

  % scale averaging time constant based on number of samples 

  tc = 0.95*Ndft/Fs;
  %tc = .95;
  % Update mag DFT  ---------------------------------------------

  numffts = floor(length(sf)/Ndft);
  h = hanning(Ndft);
  for i=1:numffts
    a = (i-1)*Ndft+1; b = i*Ndft;
    Sf = abs(fft(sf(a:b) .* h, Ndft));
    Sf(1:st) = 0; Sf(en:Ndft/2) = 0;
    states.Sf = (1-tc)*states.Sf + tc*Sf(1:Ndft/2);
  end

  f = []; a = [];
  Sf = states.Sf;

  %figure(8)
  %clf
  %plot(Sf(1:Ndft/2));

  % Search for each tone --------------------------------------------------------

  for m=1:ntones
    [tone_amp tone_index] = max(Sf(1:Ndft/2));

    f = [f (tone_index-1)*Fs/Ndft];
    a = [a tone_amp];

    % zero out region min_tone_spacing/2 either side of max so we can find next highest peak
    % closest spacing for non-coh mFSK is Rs

    st = tone_index - floor((min_tone_spacing/2)*Ndft/Fs);
    st = max(1,st);
    en = tone_index + floor((min_tone_spacing/2)*Ndft/Fs);
    en = min(Ndft/2,en);
    Sf(st:en) = 0;
  end

  states.f = sort(f);
end


% Given a buffer of nin input Rs baud FSK samples, returns nsym bits.
%
% nin is the number of input samples required by demodulator.  This is
% time varying.  It will nominally be N (8000), and occasionally N +/- 
% Ts/2 (e.g. 8080 or 7920).  This is how we compensate for differences between the
% remote tx sample clock and our sample clock.  This function always returns
% N/Ts (50) demodulated bits.  Variable number of input samples, constant number
% of output bits.

function [rx_bits states] = fsk_horus_demod(states, sf)
  M = states.M;
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
  f = states.f;

  assert(length(sf) == nin);

  % down convert and filter at rate P ------------------------------

  % update filter (integrator) memory by shifting in nin samples
  
  nold = Nmem-nin; % number of old samples we retain

  f_dc = states.f_dc; 
  f_dc(:,1:nold) = f_dc(:,Nmem-nold+1:Nmem);

  % shift down to around DC, ensuring continuous phase from last frame

  for m=1:M
    phi_vec = states.phi(m) + (1:nin)*2*pi*f(m)/Fs;
    f_dc(m,nold+1:Nmem) = sf .* exp(j*phi_vec)';
    states.phi(m)  = phi_vec(nin);
    states.phi(m) -= 2*pi*floor(states.phi(m)/(2*pi));
  end

  % save filter (integrator) memory for next time

  states.f_dc = f_dc;

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
    for m=1:M
      f_int(m,i) = sum(f_dc(m,st:en));
    end
  end
  states.f_int = f_int;

  % fine timing estimation -----------------------------------------------

  % Non linearity has a spectral line at Rs, with a phase
  % related to the fine timing offset.  See:
  %   http://www.rowetel.com/blog/?p=3573 
  % We have sampled the integrator output at Fs=P samples/symbol, so
  % lets do a single point DFT at w = 2*pi*f/Fs = 2*pi*Rs/(P*Rs)
  %
  % Note timing non-lineariry derived by experiment.  Not quite sure what I'm doing here.....
  % but it gives 0dB impl loss for 2FSK Eb/No=9dB, testmode 1:
  %   Fs: 8000 Rs: 50 Ts: 160 nsym: 50
  %   frames: 200 Tbits: 9700 Terrs: 93 BER 0.010

  Np = length(f_int(1,:));
  w = 2*pi*(Rs)/(P*Rs);
  timing_nl = sum(abs(f_int(:,:)).^2);
  x = timing_nl * exp(-j*w*(0:Np-1))';
  norm_rx_timing = angle(x)/(2*pi);
  %norm_rx_timing = 0;
  rx_timing = norm_rx_timing*P;

  states.x = x;
  states.timing_nl = timing_nl;
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

  f_int_resample = zeros(M,nsym);
  rx_bits = zeros(1,nsym);
  tone_max = rx_bits_sd = zeros(1,nsym);

  for i=1:nsym
    st = i*P+1;
    f_int_resample(:,i) = f_int(:,st+low_sample)*(1-fract) + f_int(:,st+high_sample)*fract;
    %rx_bits(i) = abs(f_int_resample(2,i)) > abs(f_int_resample(1,i));
    %rx_bits_sd(i) = abs(f_int_resample(2,i)) - abs(f_int_resample(2,i));
    [tone_max(i) tone_index] = max(f_int_resample(:,i));
    if M == 2
      rx_bits(i) = tone_index - 1;
      rx_bits_sd(i) = abs(f_int_resample(2,i)) - abs(f_int_resample(2,i));
    else
      demap = [[0 0]; [0 1]; [1 0]; [1 1]];      
      rx_bits(2*i-1:2*i) = demap(tone_index,:);
    end

  end

  states.f_int_resample = f_int_resample;
  states.rx_bits_sd = rx_bits_sd;

  % Eb/No estimation

  tone_max = abs(tone_max);
  states.EbNodB = -6 + 20*log10(1E-6+mean(tone_max)/(1E-6+std(tone_max)));
endfunction


% Look for unique word and return index of first UW bit, or -1 if no
% UW found Sometimes there may be several matches, returns the
% position of the best match to UW.

function [uw_start best_corr corr] = find_uw(states, start_bit, rx_bits)
  uw = states.uw;

  mapped_rx_bits = 2*rx_bits - 1;
  best_corr = 0;
  uw_start = -1;
  found_uw = 0;

  % first first UW in buffer that exceeds threshold
  
  for i=start_bit:length(rx_bits) - length(uw)
    corr(i)  = mapped_rx_bits(i:i+length(uw)-1) * uw';
    if (found_uw == 0) && (corr(i) >= states.uw_thresh)
      uw_start = i;
      best_corr = corr;
      found_uw = 1;
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

  uw_loc = find_uw(states.rtty, bit, rx_bits_log, states.verbose);
  
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
      end

      % update memory of previous packet, we use this to guess where errors may be
      if crc_ok || crc_flipped_ok
        states.prev_pkt = rx_bits_log(uw_loc+length(states.rtty.uw):uw_loc+states.rtty.max_packet_len);
      end

      if crc_ok
        str = sprintf("%s CRC OK", str);
      else
        if crc_flipped_ok
          str = sprintf("%s fixed", str_flipped);
        else
          str = sprintf("%s CRC BAD", str);
        end
      end
      printf("%s\n", str);
    end

    % look for next packet

    bit = uw_loc + length(states.rtty.uw);
    uw_loc = find_uw(states.rtty, bit, rx_bits_log, states.verbose);

  endwhile
endfunction
 

% Extract as many binary packets as we can from a great big buffer of bits,
% and send them to the C decoder for FEC decoding.
% horus_l2 can be compiled a bunch of different ways.  You need to
% compile with:
%   codec2-dev/src$ gcc horus_l2.c -o horus_l2 -Wall -DDEC_RX_BITS -DHORUS_L2_RX

function corr_log = extract_and_decode_binary_packets(states, rx_bits_log)
  corr_log = [];

  % use UWs to delimit start and end of data packets

  bit = 1;
  nbits = length(rx_bits_log);

  [uw_loc best_corr corr] = find_uw(states.binary, bit, rx_bits_log, states.verbose);
  corr_log = [corr_log corr];

  while (uw_loc != -1)

    if (uw_loc+states.binary.max_packet_len) < nbits
      % printf("uw_loc: %d best_corr: %d\n", uw_loc, best_corr);

      % OK we have a packet delimited by two UWs.  Lets convert the bit
      % stream into bytes and save for decoding

      pin = uw_loc;    
      for i=1:45
        rx_bytes(i) = rx_bits_log(pin:pin+7) * (2.^(7:-1:0))';
        pin += 8;
        %printf("%d 0x%02x\n", i, rx_bytes(i));
      end

      f=fopen("horus_rx_bits_binary.bin","wb");
      fwrite(f, rx_bytes, "uchar");
      fclose(f);

      % optionally write packet to disk to use as horus_tx_bits_binary.txt
      f=fopen("horus_rx_bits_binary.txt","wt");
      for i=uw_loc:uw_loc+45*8-1
        fprintf(f, "%d ", rx_bits_log(i));
      end
      fclose(f);

      system("../src/horus_l2");  % compile instructions above
    end

    bit = uw_loc + length(states.binary.uw);
    [uw_loc best_corr corr] = find_uw(states.binary, bit, rx_bits_log, states.verbose);
    corr_log = [corr_log corr];
   
  endwhile
endfunction
 

% BER counter and test frame sync logic

function states = ber_counter(states, test_frame, rx_bits_buf)
  nbit = states.nbit;
  state = states.ber_state;
  next_state = state;

  if state == 0

    % try to sync up with test frame

    nerrs_min = nbit;
    for i=1:nbit
      error_positions = xor(rx_bits_buf(i:nbit+i-1), test_frame);
      nerrs = sum(error_positions);
      if nerrs < nerrs_min
        nerrs_min = nerrs;
        states.coarse_offset = i;
      end
    end
    if nerrs_min/nbit < 0.05 
      next_state = 1;
    end
    if bitand(states.verbose,0x4)
      printf("coarse offset: %d nerrs_min: %d next_state: %d\n", states.coarse_offset, nerrs_min, next_state);
    end
  end

  if state == 1  

    % we're synced up, lets measure bit errors

    error_positions = xor(rx_bits_buf(states.coarse_offset:states.coarse_offset+nbit-1), test_frame);
    nerrs = sum(error_positions);
    if nerrs/nbit > 0.1
      next_state = 0;
    else
      states.Terrs += nerrs;
      states.Tbits += nbit;
      states.nerr_log = [states.nerr_log nerrs];
    end
  end

  states.ber_state = next_state;
endfunction


% Alternative stateless BER counter that works on packets that may have gaps betwene them

function states = ber_counter_packet(states, test_frame, rx_bits_buf)
  ntestframebits = states.ntestframebits;
  nbit = states.nbit;

  % look for offset with min errors

  nerrs_min = ntestframebits; coarse_offset = 1;
  for i=1:nbit
    error_positions = xor(rx_bits_buf(i:ntestframebits+i-1), test_frame);
    nerrs = sum(error_positions);
    %printf("i: %d nerrs: %d\n", i, nerrs);
    if nerrs < nerrs_min
      nerrs_min = nerrs;
      coarse_offset = i;
    end
  end

  % if less than threshold count errors

  if nerrs_min/ntestframebits < 0.05 
    states.Terrs += nerrs_min;
    states.Tbits += ntestframebits;
    states.nerr_log = [states.nerr_log nerrs_min];
    if bitand(states.verbose, 0x4)
      printf("coarse_offset: %d nerrs_min: %d\n", coarse_offset, nerrs_min);
    end
  end
endfunction


% simulation of tx and rx side, add noise, channel impairments ----------------------
%
% test_frame_mode     Description
% 1                   BER testing using known test frames
% 2                   random bits
% 3                   repeating sequence of all symbols
% 4                   Horus RTTY
% 5                   Horus Binary
% 6                   Horus High Speed: A 8x oversampled modem, e.g. Fs=9600, Rs=1200
%                     which is the same as Fs=921600 Rs=115200
%                     Uses packed based BER counter

function run_sim(test_frame_mode, frames = 10, EbNodB = 100)
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

  if test_frame_mode < 4
    % horus rtty config ---------------------
    states = fsk_horus_init(8000, 50, 4);
    %states = fsk_horus_init_hbr(8000, 10, 400, 4); % EME
  end

  if test_frame_mode == 4
    % horus rtty config ---------------------
    states = fsk_horus_init(8000, 100);
    states.tx_bits_file = "horus_tx_bits_rtty.txt"; % Octave file of bits we FSK modulate
  end
                               
  if test_frame_mode == 5
    % horus binary config ---------------------
    states = fsk_horus_init(8000, 50, 4);
    states.tx_bits_file = "horus_tx_bits_binary.txt"; % Octave file of bits we FSK modulate
  end

  if test_frame_mode == 6
    % horus high speed ---------------------
    states = fsk_horus_init_hbr(9600, 8, 1200, 2, 16);
    states.tx_bits_file = "horus_high_speed.bin";
  end

  % Tones must be at least Rs apart for ideal non-coherent FSK

  if states.M == 2
    states.ftx = 1200 + [ 0 states.Rs ];
  else
    states.ftx = 1200 + 2*states.Rs*(1:4)
    %states.ftx = 200 + states.Rs*(1:4); % EME
  end

  % ----------------------------------------------------------------------

  states.verbose = 0x1;
  M = states.M;
  N = states.N;
  P = states.P;
  Rs = states.Rs;
  nsym = states.nsym;
  nbit = states.nbit;
  Fs = states.Fs;
  states.df(1:M) = df;
  states.dA(1:M) = dA;

  % optional noise.  Useful for testing performance of waveforms from real world modulators

  EbNo = 10^(EbNodB/10);
  variance = states.Fs/(states.Rs*EbNo*states.bitspersymbol);

  % set up tx signal with payload bits based on test mode

  if (test_frame_mode == 1)
    % test frame of bits, which we repeat for convenience when BER testing
    states.ntestframebits = states.nbit;
    test_frame = round(rand(1, states.ntestframebits));
    tx_bits = [];
    for i=1:frames+1
      tx_bits = [tx_bits test_frame];
    end
  end

  if test_frame_mode == 2
    % random bits, just to make sure sync algs work on random data
    tx_bits = round(rand(1, states.nbit*(frames+1)));
  end

  if test_frame_mode == 3
    % repeating sequence of all symbols
    % great for initial test of demod if nothing else works, 
    % look for this pattern in rx_bits
    if M == 2
       % ...10101...
      tx_bits = zeros(1, states.nbit*(frames+1));
      tx_bits(1:2:length(tx_bits)) = 1;
    else
      % repeat each possible 4fsk symbol
      pattern = [0 0 0 1 1 0 1 1];
      %pattern = [0 0 0 1 1 1 1 0];
      nrepeats = states.nbit*(frames+1)/length(pattern);
      tx_bits = [];
      for b=1:nrepeats
        tx_bits = [tx_bits pattern];
      end   
      %tx_bits = zeros(1, states.nbit*(frames+1));
    end
  end
 
  if (test_frame_mode == 4) || (test_frame_mode == 5)

    % load up a horus msg from disk and modulate that

    test_frame = load(states.tx_bits_file);
    ltf = length(test_frame);
    ntest_frames = ceil((frames+1)*nbit/ltf);
    tx_bits = [];
    for i=1:ntest_frames
      tx_bits = [tx_bits test_frame];
    end
  end

  if test_frame_mode == 6
    states.verbose += 0x4;
    ftmp = fopen(states.tx_bits_file, "rb"); test_frame = fread(ftmp,Inf,"char")'; fclose(ftmp);
    states.ntestframebits = length(test_frame);
    printf("length test frame: %d\n", states.ntestframebits);
    %test_frame = rand(1,states.ntestframebits) > 0.5;

    tx_bits = [];
    for i=1:frames+1
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
  printf("SNRdB meas: %4.1f\n", 10*log10(var(tx)/var(noise)));

  % dump simulated rx file

  ftx=fopen("fsk_horus.raw","wb"); rxg = rx*1000; fwrite(ftx, rxg, "short"); fclose(ftx);

  timing_offset_samples = round(timing_offset*states.Ts);
  st = 1 + timing_offset_samples;
  rx_bits_buf = zeros(1,nbit+states.ntestframebits);
  x_log = [];
  timing_nl_log = [];
  norm_rx_timing_log = [];
  f_int_resample_log = [];
  f_log = [];
  EbNodB_log = [];
  rx_bits_log = [];
  rx_bits_sd_log = [];

  % main loop ---------------------------------------------------------------

  run_frames = floor(length(rx)/N)-1;
  for f=1:run_frames

    % extract nin samples from input stream

    nin = states.nin;
    en = st + states.nin - 1;

    if en < length(rx) % due to nin variations its possible to overrun buffer
      sf = rx(st:en);
      st += nin;

      % demodulate to stream of bits

      states = est_freq(states, sf, states.M);
      %states.f = [1200 1400 1600 1800];
      [rx_bits states] = fsk_horus_demod(states, sf);

      rx_bits_buf(1:states.ntestframebits) = rx_bits_buf(nbit+1:states.ntestframebits+nbit);
      rx_bits_buf(states.ntestframebits+1:states.ntestframebits+nbit) = rx_bits;
      %rx_bits_buf(1:nbit) = rx_bits_buf(nbit+1:2*nbit);
      %rx_bits_buf(nbit+1:2*nbit) = rx_bits;

      rx_bits_log = [rx_bits_log rx_bits];
      rx_bits_sd_log = [rx_bits_sd_log states.rx_bits_sd];

      norm_rx_timing_log = [norm_rx_timing_log states.norm_rx_timing];
      x_log = [x_log states.x];
      timing_nl_log = [timing_nl_log states.timing_nl];
      f_int_resample_log = [f_int_resample_log abs(states.f_int_resample(:,:))];
      f_log = [f_log; states.f];
      EbNodB_log = [EbNodB_log states.EbNodB];

      if test_frame_mode == 1
        states = ber_counter(states, test_frame, rx_bits_buf);
      end
      if test_frame_mode == 6
        states = ber_counter_packet(states, test_frame, rx_bits_buf);
      end
   end
  end

  % print stats, count errors, decode packets  ------------------------------------------

  if (test_frame_mode == 1) || (test_frame_mode == 6)
    printf("frames: %d EbNo: %3.2f Tbits: %d Terrs: %d BER %4.3f\n", frames, EbNodB, states.Tbits,states. Terrs, states.Terrs/states.Tbits);
  end

  if test_frame_mode == 4
    extract_and_print_rtty_packets(states, rx_bits_log, rx_bits_sd_log)
  end

  if test_frame_mode == 5
    extract_and_decode_binary_packets(states, rx_bits_log);
  end

  figure(1);
  plot(f_int_resample_log','+')
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
  axis([1 run_frames -1 1])
  title('norm fine timing')
  subplot(212)
  plot(states.nerr_log)
  title('num bit errors each frame')

  figure(4)
  clf
  subplot(211)
  one_sec_rx = rx(1:min(Fs,length(rx)));
  plot(one_sec_rx)
  title('rx signal at demod input')
  subplot(212)
  plot(abs(fft(one_sec_rx)))

  figure(5)
  clf
  plot(f_log,'+')
  title('tone frequencies')
  axis([1 run_frames 0 Fs/2])

  figure(6)
  clf
  plot(EbNodB_log);
  title('Eb/No estimate')

  figure(7)
  clf
  subplot(211)
  X = abs(fft(timing_nl_log));
  plot(X(1:length(X)/2))
  subplot(212)
  plot(abs(timing_nl_log(1:100)))

 endfunction


% demodulate a file of 8kHz 16bit short samples --------------------------------

function rx_bits_log = demod_file(filename, test_frame_mode, noplot=0, EbNodB=100)
  fin = fopen(filename,"rb"); 
  more off;

  %states = fsk_horus_init(96000, 1200);

  if test_frame_mode == 4
    % horus rtty config ---------------------
    states = fsk_horus_init(8000, 100, 2);
    uwstates = fsk_horus_init_rtty_uw(states);
    states.ntestframebits = states.nbits;
  end
                               
  if test_frame_mode == 5
    % horus binary config ---------------------
    states = fsk_horus_init(8000, 50, 4);
    uwstates = fsk_horus_init_binary_uw;
    states.ntestframebits = states.nbits;
  end

  states.verbose = 0x1 + 0x8;

  if test_frame_mode == 6
    % Horus high speed config --------------
    states = fsk_horus_init_hbr(9600, 8, 1200, 2, 16);
    states.tx_bits_file = "horus_high_speed.bin";
    states.verbose += 0x4;
    ftmp = fopen(states.tx_bits_file, "rb"); test_frame = fread(ftmp,Inf,"char")'; fclose(ftmp);
    states.ntestframebits = length(test_frame);
    printf("length test frame: %d\n", states.ntestframebits);
  end

  if test_frame_mode == 7
    % 800XA 4FSK modem --------------
    states = fsk_horus_init_hbr(8000, 10, 400, 4, 256);
    states.tx_bits_file = "horus_high_speed.bin";
    states.verbose += 0x4;
    ftmp = fopen(states.tx_bits_file, "rb"); test_frame = fread(ftmp,Inf,"char")'; fclose(ftmp);
    states.ntestframebits = length(test_frame);
    printf("length test frame: %d\n", states.ntestframebits);
  end

  N = states.N;
  P = states.P;
  Rs = states.Rs;
  nsym = states.nsym;
  nbit = states.nbit;

  frames = 0;
  rx = [];
  rx_bits_log = [];
  rx_bits_sd_log = [];
  norm_rx_timing_log = [];
  f_int_resample_log = [];
  EbNodB_log = [];
  ppm_log = [];
  f_log = [];
  rx_bits_buf = zeros(1,nbit + states.ntestframebits);

  % optional noise.  Useful for testing performance of waveforms from real world modulators

  EbNo = 10^(EbNodB/10);
  ftmp = fopen(filename,"rb"); s = fread(ftmp,Inf,"short"); fclose(ftmp); tx_pwr = var(s);
  variance = (tx_pwr/2)*states.Fs/(states.Rs*EbNo*states.bitspersymbol);

  % First extract raw bits from samples ------------------------------------------------------

  printf("demod of raw bits....\n");

  finished = 0;
  while (finished == 0)

    % extract nin samples from input stream

    nin = states.nin;
    [sf count] = fread(fin, nin, "short");
    rx = [rx; sf];
    
    % add optional noise

    if count
      noise = sqrt(variance)*randn(count,1);
      sf += noise;
    end

    if count == nin
      frames++;

      % demodulate to stream of bits

      states = est_freq(states, sf, states.M);
      %states.f = [1450 1590 1710 1850];
      [rx_bits states] = fsk_horus_demod(states, sf);

      rx_bits_buf(1:states.ntestframebits) = rx_bits_buf(nbit+1:states.ntestframebits+nbit);
      rx_bits_buf(states.ntestframebits+1:states.ntestframebits+nbit) = rx_bits;

      rx_bits_log = [rx_bits_log rx_bits];
      rx_bits_sd_log = [rx_bits_sd_log states.rx_bits_sd];
      norm_rx_timing_log = [norm_rx_timing_log states.norm_rx_timing];
      f_int_resample_log = [f_int_resample_log abs(states.f_int_resample)];
      EbNodB_log = [EbNodB_log states.EbNodB];
      ppm_log = [ppm_log states.ppm];
      f_log = [f_log; states.f];

      if test_frame_mode == 1
        states = ber_counter(states, test_frame, rx_bits_buf);
        if states.ber_state == 1
          states.verbose = 0;
        end
      end
      if test_frame_mode == 6
        states = ber_counter_packet(states, test_frame, rx_bits_buf);
      end
     else      
      finished = 1;
    end
  end
  fclose(fin);

  if noplot == 0
    printf("plotting...\n");

    figure(1);
    plot(f_log);
    hold off;

    figure(2);
    plot(f_int_resample_log','+')

    figure(3)
    clf
    subplot(211)
    plot(norm_rx_timing_log)
    axis([1 frames -0.5 0.5])
    title('norm fine timing')
    grid
    subplot(212)
    plot(states.nerr_log)
    title('num bit errors each frame')
 
    figure(4)
    clf
    plot(EbNodB_log);
    title('Eb/No estimate')

    figure(5)
    clf
    rx_nowave = rx(1000:length(rx));
    subplot(211)
    plot(rx_nowave(1:states.Fs));
    title('input signal to demod (1 sec)')
    xlabel('Time (samples)');
    axis([1 states.Fs -35000 35000])

    % normalise spectrum to 0dB full scale with a 32767 sine wave input

    subplot(212)
    RxdBFS = 20*log10(abs(fft(rx_nowave(1:states.Fs)))) - 20*log10((states.Fs/2)*32767);
    plot(RxdBFS)
    axis([1 states.Fs/2 -80 0])
    xlabel('Frequency (Hz)');

    figure(6);
    clf
    plot(ppm_log)
    title('Sample clock (baud rate) offset in PPM');
  end

  if (test_frame_mode == 1) || (test_frame_mode == 6)
    printf("frames: %d Tbits: %d Terrs: %d BER %4.3f EbNo: %3.2f\n", frames, states.Tbits,states. Terrs, states.Terrs/states.Tbits, mean(EbNodB_log));
  end

  % we can decode both protocols at the same time

  if (test_frame_mode == 4) || (test_frame_mode == 5)
    extract_and_print_rtty_packets(states, rx_bits_log, rx_bits_sd_log)
    corr_log = extract_and_decode_binary_packets(states, rx_bits_log);

    figure(8);
    clf
    plot(corr_log);
    hold on;
    plot([1 length(corr_log)],[states.binary.uw_thresh states.binary.uw_thresh],'g');
    hold off;
    title('UW correlation');
  end

endfunction


% run test functions from here during development

if exist("fsk_horus_as_a_lib") == 0
  %run_sim(1, 100, 6);
  %rx_bits = demod_file("~/Desktop/115.wav",6,0,90);
  %rx_bits = demod_file("fsk_horus.raw",5);
  %rx_bits = demod_file("~/Desktop/4FSK_Binary_NoLock.wav",4);
  %rx_bits = demod_file("~/Desktop/phorus_binary_ascii.wav",4);
  %rx_bits = demod_file("~/Desktop/binary/horus_160102_binary_rtty_2.wav",4);
  %rx_bits = demod_file("~/Desktop/horus_160102_vk5ei_capture2.wav",4);
  %rx_bits = demod_file("~/Desktop/horus_rtty_binary.wav",4);
  %rx_bits = demod_file("~/Desktop/FSK_4FSK.wav",4);
  %rx_bits = demod_file("t.raw",5);
  %rx_bits = demod_file("~/Desktop/fsk_horus_10dB_1000ppm.wav",4);
  %rx_bits = demod_file("~/Desktop/fsk_horus_6dB_0ppm.wav",4);
  %rx_bits = demod_file("test.raw",1,1);
  %rx_bits = demod_file("/dev/ttyACM0",1);
  %rx_bits = demod_file("fsk_horus_rx_1200_96k.raw",1);
  %rx_bits = demod_file("mp.raw",4);
  %rx_bits = demod_file("~/Desktop/launchbox_v2_landing_8KHz_final.wav",4);
  rx_bits = demod_file("~/Desktop/bench_test_002.wav",7);
end
