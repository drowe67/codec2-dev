% ofdm_rx.m
% David Rowe April 2017
%
% OFDM file based rx.

#{
  TODO: 
    [ ] proper EsNo estimation
    [ ] some sort of real time GUI display to watch signal evolving
    [ ] est SNR or Eb/No of recieved signal
    [ ] way to fall out of sync
#}

function ofdm_rx(filename, interleave_frames = 1)
  ofdm_lib;
  ldpc;
  more off;

  % init modem

  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 16; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;
  states.verbose = 1;

  % Set up LDPC code

  mod_order = 4; bps = 2; modulation = 'QPSK'; mapping = 'gray';
  demod_type = 0; decoder_type = 0; max_iterations = 100;

  EsNo = 10; % TODO: fixme

  init_cml('/home/david/Desktop/cml/');
  load HRA_112_112.txt
  [code_param framesize rate] = ldpc_init_user(HRA_112_112, modulation, mod_order, mapping);
  assert(Nbitsperframe == code_param.code_bits_per_frame);

  % load real samples from file

  Ascale= 2E5;
  frx=fopen(filename,"rb"); rx = 2*fread(frx, Inf, "short")/4E5; fclose(frx);
  Nsam = length(rx); Nframes = floor(Nsam/Nsamperframe);
  prx = 1;

  % buffers for interleaved frames

  Nsymbolsperframe = code_param.code_bits_per_frame/bps
  Nsymbolsperinterleavedframe = interleave_frames*Nsymbolsperframe
  rx_np = zeros(1, Nsymbolsperinterleavedframe);
  rx_amp = zeros(1, Nsymbolsperinterleavedframe);

  % OK generate tx frame for BER calcs

  rand('seed', 100);
  tx_bits = []; tx_codewords = [];
  for f=1:interleave_frames
    atx_bits = round(rand(1,code_param.data_bits_per_frame));
    tx_bits = [tx_bits atx_bits];
    acodeword = LdpcEncode(atx_bits, code_param.H_rows, code_param.P_matrix);
    tx_codewords = [tx_codewords acodeword];
  end

  % used for rx frame sync on interleaved symbols - we demod the
  % entire interleaved frame to raw bits

  tx_symbols = [];
  for s=1:Nsymbolsperinterleavedframe
    tx_symbols = [tx_symbols qpsk_mod( tx_codewords(2*(s-1)+1:2*s) )];
  end

  tx_symbols = gp_interleave(tx_symbols);

  tx_bits_raw = [];
  for s=1:Nsymbolsperinterleavedframe
    tx_bits_raw = [tx_bits_raw qpsk_demod(tx_symbols(s))];
  end

  % init logs and BER stats

  rx_bits = []; rx_np_log = []; timing_est_log = []; delta_t_log = []; foff_est_hz_log = [];
  phase_est_pilot_log = [];
  Terrs = Tbits = Terrs_coded = Tbits_coded = Tpackets = Tpacketerrs = 0;
  Nbitspervocframe = 28;

  % 'prime' rx buf to get correct coarse timing (for now)

  prx = 1;
  nin = Nsamperframe+2*(M+Ncp);
  states.rxbuf(Nrxbuf-nin+1:Nrxbuf) = rx(prx:nin);
  prx += nin;

  state = 'searching'; frame_count = 0;

  % main loop ----------------------------------------------------------------

  for f=1:Nframes

    % insert samples at end of buffer, set to zero if no samples
    % available to disable phase estimation on future pilots on last
    % frame of simulation

    lnew = min(Nsam-prx,states.nin);
    rxbuf_in = zeros(1,states.nin);

    if lnew
      rxbuf_in(1:lnew) = rx(prx:prx+lnew-1);
    end
    prx += states.nin;

    [rx_bits_raw states aphase_est_pilot_log arx_np arx_amp] = ofdm_demod(states, rxbuf_in);
    frame_count++;

    % update sliding windows of rx-ed symbols and symbol amplitudes

    rx_np(1:Nsymbolsperinterleavedframe-Nsymbolsperframe) = rx_np(Nsymbolsperframe+1:Nsymbolsperinterleavedframe);
    rx_np(Nsymbolsperinterleavedframe-Nsymbolsperframe+1:Nsymbolsperinterleavedframe) = arx_np;
    rx_amp(1:Nsymbolsperinterleavedframe-Nsymbolsperframe) = rx_amp(Nsymbolsperframe+1:Nsymbolsperinterleavedframe);
    rx_amp(Nsymbolsperinterleavedframe-Nsymbolsperframe+1:Nsymbolsperinterleavedframe) = arx_amp;

    printf("f: %d state: %s frame_count: %d\n", f, state, frame_count);

    % If looking for sync: check raw BER on frame just received
    % against all possible positions in the interleaver frame.

    % iterate state machine ------------------------------------

    next_state = state;
    if strcmp(state,'searching')  

      % If looking for sync: check raw BER on frame just received
      % against all possible positions in the interleaver frame.

      for ff=1:interleave_frames
        st = (ff-1)*Nbitsperframe+1; en = st + Nbitsperframe - 1;
        errors = xor(tx_bits_raw(st:en), rx_bits_raw);
        Nerrs = sum(errors); 
        printf("  ff: %d Nerrs: %d\n", ff, Nerrs);
        if Nerrs/Nbitsperframe < 0.1
          next_state = 'synced';
          frame_count = ff;
          if f < interleave_frames
            % point trying a LDPC decode if we don't have a full frame!
            frame_count -= interleave_frames;
          end
        end
      end
    end
    state = next_state;

    if strcmp(state,'searching') 

      % still searching? Attempt coarse timing estimate (i.e. detect start of frame)

      st = M+Ncp + Nsamperframe + 1; en = st + 2*Nsamperframe; 
      [ct_est foff_est] = coarse_sync(states, states.rxbuf(st:en), states.rate_fs_pilot_samples);
      if states.verbose
        printf("  Nerrs: %d ct_est: %4d foff_est: %3.1f\n", Nerrs, ct_est, foff_est);
      end

      % calculate number of samples we need on next buffer to get into sync
     
      states.nin = Nsamperframe + ct_est - 1;

      % reset modem states

      states.sample_point = states.timing_est = 1;
      states.foff_est_hz = foff_est;
    end
    
    if strcmp(state,'synced') && (frame_count == interleave_frames)

      % de-interleave QPSK symbols

      arx_np = gp_deinterleave(rx_np);
      arx_amp = gp_deinterleave(rx_amp);

      % LDPC decode

      rx_bits = [];
      for ff=1:interleave_frames
        st = (ff-1)*Nbitsperframe/bps+1; en = st + Nbitsperframe/bps - 1;
        rx_codeword = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, arx_np(st:en), min(EsNo,30), arx_amp(st:en));
        rx_bits = [rx_bits rx_codeword(1:code_param.data_bits_per_frame)];
      end

      errors_coded = xor(tx_bits, rx_bits);
      Nerrs_coded = sum(errors_coded);
      Terrs_coded += Nerrs_coded;
      Tbits_coded += code_param.data_bits_per_frame*interleave_frames;
      Nerrs_coded_log(f) = Nerrs_coded;

      printf("  Nerrs_coded: %d\n", Nerrs_coded);

      % we are in sync so log states and bit/packet error stats

      rx_np_log = [rx_np_log arx_np];
      timing_est_log = [timing_est_log states.timing_est];
      delta_t_log = [delta_t_log states.delta_t];
      foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
      phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];

      % measure uncoded bit errors

      rx_bits_raw = [];
      for s=1:Nsymbolsperinterleavedframe
        rx_bits_raw = [rx_bits_raw qpsk_demod(rx_np(s))];
      end
      errors = xor(tx_bits_raw, rx_bits_raw);
      Nerrs = sum(errors);
      Terrs += Nerrs;
      Nerrs_log(f) = Nerrs;
      Tbits += code_param.code_bits_per_frame*interleave_frames;

      % measure packet errors based on Codec 2 packet size

      Nvocframes = floor(code_param.data_bits_per_frame*interleave_frames/Nbitspervocframe);
      for fv=1:Nvocframes
        st = (fv-1)*Nbitspervocframe + 1;
        en = st + Nbitspervocframe - 1;
        Nvocpacketerrs = sum(xor(tx_bits(st:en), rx_bits(st:en)));
        if Nvocpacketerrs
          Tpacketerrs++;
        end
        Tpackets++;
      end

      frame_count = 0;
    end
  end

  printf("Coded BER: %5.4f Tbits: %5d Terrs: %5d PER: %5.4f Tpacketerrs: %5d Tpackets: %5d\n", 
         Terrs_coded/Tbits_coded, Tbits_coded, Terrs_coded, Tpacketerrs/Tpackets, Tpacketerrs, Tpackets);
  printf("Raw BER..: %5.4f Tbits: %5d Terrs: %5d\n", Terrs/Tbits, Tbits, Terrs);

  figure(1); clf; 
  plot(rx_np_log,'+');
  axis([-2 2 -2 2]);
  title('Scatter');

  figure(2); clf;
  plot(phase_est_pilot_log(:,2:Nc+1),'g+', 'markersize', 5); 
  title('Phase est');
  axis([1 length(phase_est_pilot_log) -pi pi]);  

  figure(3); clf;
  subplot(211)
  stem(delta_t_log)
  title('delta t');
  subplot(212)
  plot(timing_est_log);
  title('timing est');

  figure(4); clf;
  plot(foff_est_hz_log)
  axis([1 max(Nframes,2) -2 2]);
  title('Fine Freq');
  ylabel('Hz')

  figure(5); clf;
  subplot(211)
  title('Nerrs Log')
  stem(Nerrs_log);
  subplot(212)
  stem(Nerrs_coded_log);
endfunction
