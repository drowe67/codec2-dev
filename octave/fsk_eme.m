% fsk_eme.m
% David Rowe July 2017
%
% FSK modem for digital voice over 1296 MHz Earth-Moon-Earth (EME) links

fsk_lib;

function run_sim(frames = 10, EbNodB = 100)
  Fs = 8000; Rs = 8; M = 16; Fsep = 50;

  timing_offset = 0.0; % see resample() for clock offset below
  fading = 0;          % modulates tx power at 2Hz with 20dB fade depth, 
                       % to simulate balloon rotating at end of mission
  df     = 0;          % tx tone freq drift in Hz/s
  dA     = 1;          % amplitude imbalance of tones (note this affects Eb so not a gd idea)

  more off
  rand('state',1); 
  randn('state',1);

  states = fsk_init(Fs, 50, M);

  states.ftx = 1000 + Fsep*(-M/2:M/2-1);

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

  % test frame of bits, which we repeat for convenience when BER testing

  states.ntestframebits = states.nbit;
  test_frame = round(rand(1, states.ntestframebits));
  tx_bits = [];
  for i=1:frames+1
    tx_bits = [tx_bits test_frame];
  end


  tx = fsk_mod(states, tx_bits);

  %tx = resample(tx, 1000, 1001); % simulated 1000ppm sample clock offset

  if fading
     ltx = length(tx);
     tx = tx .* (1.1 + cos(2*pi*2*(0:ltx-1)/Fs))'; % min amplitude 0.1, -20dB fade, max 3dB
  end

  noise = sqrt(variance)*randn(length(tx),1);
  rx    = tx + noise;
  printf("SNRdB meas: %4.1f\n", 10*log10(var(tx)/var(noise)));

  % dump simulated rx file

  g = 1000; 
  ftx=fopen("fsk_horus_tx.raw","wb"); txg = tx*g; fwrite(ftx, txg, "short"); fclose(ftx);
  frx=fopen("fsk_horus_rx.raw","wb"); rxg = rx*g; fwrite(frx, rxg, "short"); fclose(frx);

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
      states.f = states.ftx;
      [rx_bits states] = fsk_demod(states, sf);

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

      states = ber_counter(states, test_frame, rx_bits_buf);
    end
  end

  % print stats, count errors, decode packets  ------------------------------------------

  printf("frames: %d EbNo: %3.2f Tbits: %d Terrs: %d BER %4.3f\n", frames, EbNodB, states.Tbits,states. Terrs, states.Terrs/states.Tbits);

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


% Start simulations here -------------------------------------------------------------

more off; format;

%run_sim(1, 2, 100, 9);
run_sim(10, 20);
