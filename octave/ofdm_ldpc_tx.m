% ofdm_ldpc_tx.m
% David Rowe April 2017
%
% File based ofdm tx with LDPC encoding and interleaver.  Generates a
% file of ofdm samples, including optional channel simulation.

#{
  Examples:
 
  i) 4 frame interleaver, 10 seconds, AWGN channel at (coded) Eb/No=3dB

    octave:4> ofdm_ldpc_tx('awgn_ebno_3dB_700d.raw', 4, 10,3);

  ii) 4 frame interleaver, 10 seconds, HF channel at (coded) Eb/No=6dB

    ofdm_ldpc_tx('hf_ebno_6dB_700d.raw', 4, 10, 6, 'hf');
#}


#{
  TODO: 
    [ ] measure and report raw and coded BER
    [ ] maybe 10s worth of frames, sync up to any one automatically
        + or start with maybe 10 frames
        + measure BER match on each one
    [ ] model clipping/PA compression
    [ ] sample clock offsets
    [ ] compare with same SNR from pathsim
    [ ] How to pack arbitrary frames into ofdm frame and codec 2 frames
        + integer number of ofdm frames?
        + how to sync at other end
 
#}

function ofdm_ldpc_tx(filename, interleave_frames = 1, Nsec, EbNodB=100, channel='awgn', freq_offset_Hz=0)
  ofdm_lib;
  ldpc;
  gp_interleaver;

  % init modem

  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 17; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;

  % Set up LDPC code

  mod_order = 4; bps = 2; modulation = 'QPSK'; mapping = 'gray';

  init_cml('~/cml/');
  load HRA_112_112.txt
  [code_param framesize rate] = ldpc_init_user(HRA_112_112, modulation, mod_order, mapping);
  assert(Nbitsperframe == (code_param.code_bits_per_frame + Nuwbits + Ntxtbits));

  % Generate fixed test frame of tx bits and run OFDM modulator

  Nrows = Nsec*Rs;
  Nframes = floor((Nrows-1)/Ns);

  % Adjust Nframes so we have an integer number of interleaver frames
  % in simulation

  Nframes = interleave_frames*round(Nframes/interleave_frames);

  % OK generate an interleaver frame of codewords from random data bits

  % use single frame of test bits with seed of 1
  % to keep compatable with ofdm_test_bits.c and friends
  % as per create_ldpc_test_frame
  
  rand('seed', 1);
  
  tx_bits = tx_symbols = [];
  atx_bits = round(rand(1,code_param.data_bits_per_frame));
  for f=1:interleave_frames
    tx_bits = [tx_bits atx_bits];
    codeword = LdpcEncode(atx_bits, code_param.H_rows, code_param.P_matrix);
    for b=1:2:code_param.code_bits_per_frame
      tx_symbols = [tx_symbols qpsk_mod(codeword(b:b+1))];
    end
  end
  tx_symbols = gp_interleave(tx_symbols);
  
  % generate UW and txt symbols to prepend to every frame after LDPC encoding and interleaving
 
  txt_bits = zeros(1,Ntxtbits);
  txt_symbols = [];
  for b=1:2:length(txt_bits)
    txt_symbols = [txt_symbols qpsk_mod(txt_bits(b:b+1))];
  end
  
  atx = [];
  for f=1:interleave_frames
    st = (f-1)*code_param.code_bits_per_frame/bps+1; en = st + code_param.code_bits_per_frame/bps-1;
    modem_frame = assemble_modem_frame_symbols(states, tx_symbols(st:en), txt_symbols);
    atx = [atx ofdm_txframe(states, modem_frame) ];
  end
  
  tx = [];
  for f=1:Nframes/interleave_frames
    tx = [tx atx];
  end

  % not very pretty way to process analog signals with exactly the same channel
  % todo: work out a cleaner way

  analog_hack = 0; rx_filter = 0;
  if analog_hack || rx_filter

    % simulated SSB tx filter

    [b, a] = cheby1(4, 3, [600, 2600]/(Fs/2));
    h = freqz(b,a,(600:2600)/(Fs/(2*pi)));
    filt_gain = (2600-600)/sum(abs(h) .^ 2);   % ensures power after filter == before filter
  end

  if analog_hack
    % load analog signal and convert to complex

    s = load_raw('../raw/ve9qrp_10s.raw')';
    tx = hilbert(s);

    % ssb tx filter

    tx = filter(b,a,sqrt(filt_gain)*tx);

    % normalise power to same as ofdm tx

    nom_tx_pwr = 2/(Ns*(M*M)) + Nc/(M*M);
    tx_pwr = var(tx);
    tx *= sqrt(nom_tx_pwr/tx_pwr);
  end

  Nsam = length(tx);

  % channel simulation

  EsNo = rate * bps * (10 .^ (EbNodB/10));
  variance = 1/(M*EsNo/2);
  woffset = 2*pi*freq_offset_Hz/Fs;

  SNRdB = EbNodB + 10*log10(Nc*bps*Rs*rate/3000);
  printf("EbNo: %3.1f dB  SNR(3k) est: %3.1f dB  foff: %3.1fHz inter_frms: %d ",
         EbNodB, SNRdB, freq_offset_Hz, interleave_frames);

  % set up HF model ---------------------------------------------------------------

  if strcmp(channel, 'hf')
    randn('seed',1);

    % some typical values, or replace with user supplied

    dopplerSpreadHz = 1; path_delay_ms = 1;

    path_delay_samples = path_delay_ms*Fs/1000;
    printf("Doppler Spread: %3.2f Hz Path Delay: %3.2f ms %d samples\n", dopplerSpreadHz, path_delay_ms, path_delay_samples);

    % generate same fading pattern for every run

    randn('seed',1);

    spread1 = doppler_spread(dopplerSpreadHz, Fs, (Nsec*(M+Ncp)/M)*Fs*1.1);
    spread2 = doppler_spread(dopplerSpreadHz, Fs, (Nsec*(M+Ncp)/M)*Fs*1.1);

    % sometimes doppler_spread() doesn't return exactly the number of samples we need
 
    assert(length(spread1) >= Nsam, "not enough doppler spreading samples");
    assert(length(spread2) >= Nsam, "not enough doppler spreading samples");
  end

  rx = tx;

  if strcmp(channel, 'hf')
    rx  = tx(1:Nsam) .* spread1(1:Nsam);
    rx += [zeros(1,path_delay_samples) tx(1:Nsam-path_delay_samples)] .* spread2(1:Nsam);

    % normalise rx power to same as tx

    nom_rx_pwr = 2/(Ns*(M*M)) + Nc/(M*M);
    rx_pwr = var(rx);
    rx *= sqrt(nom_rx_pwr/rx_pwr);
  end

  rx = rx .* exp(j*woffset*(1:Nsam));

  % note variance/2 as we are using real() operator, mumble,
  % reflection of -ve freq to +ve, mumble, hand wave

  noise = sqrt(variance/2)*0.5*randn(1,Nsam);
  rx = real(rx) + noise;
  printf("measured SNR: %3.2f dB\n", 10*log10(var(real(tx))/var(noise)) + 10*log10(4000) - 10*log10(3000));

  if rx_filter
    % ssb rx filter
    rx = filter(b,a,sqrt(filt_gain)*rx);
  end

  % adjusted by experiment to match rms power of early test signals

  frx=fopen(filename,"wb"); fwrite(frx, states.amp_scale*rx, "short"); fclose(frx);
endfunction
