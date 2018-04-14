% ofdm_ldpc_rx.m
% David Rowe April 2017
%
% OFDM file based rx, with LDPC and interleaver

#{
  TODO: 
    [ ] proper EsNo estimation
    [ ] some sort of real time GUI display to watch signal evolving
    [ ] est SNR or Eb/No of recieved signal
    [ ] way to fall out of sync
#}

function ofdm_ldpc_rx(filename, interleave_frames = 1, error_pattern_filename)
  ofdm_lib;
  ldpc;
  gp_interleaver;
  more off;

  % init modem

  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 17; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;
  states.verbose = 1;

  % Set up LDPC code

  mod_order = 4; bps = 2; modulation = 'QPSK'; mapping = 'gray';
  demod_type = 0; decoder_type = 0; max_iterations = 100;

  EsNo = 10; % TODO: fixme
  printf("EsNo fixed at %f - need to est from channel\n", EsNo);
  
  init_cml('/home/david/Desktop/cml/');
  load HRA_112_112.txt
  [code_param framesize rate] = ldpc_init_user(HRA_112_112, modulation, mod_order, mapping);
  assert(Nbitsperframe == (code_param.code_bits_per_frame + states.Nuwbits + states.Ntxtbits));

  % load real samples from file

  Ascale= 2E5*1.1491;
  frx=fopen(filename,"rb"); rx = 2*fread(frx, Inf, "short")/4E5; fclose(frx);
  Nsam = length(rx); Nframes = floor(Nsam/Nsamperframe);
  prx = 1;

  % buffers for interleaved frames

  Ncodedbitsperframe = code_param.code_bits_per_frame;
  Nsymbolsperframe = code_param.code_bits_per_frame/bps;
  Nuwtxtsymbolsperframe = (Nuwbits+Ntxtbits)/bps;
  Nsymbolsperinterleavedframe = interleave_frames*Nsymbolsperframe;
  rx_np = zeros(1, Nsymbolsperinterleavedframe);
  rx_amp = zeros(1, Nsymbolsperinterleavedframe);
  rx_uw = [];
  
  % OK generate tx frame for BER calcs

  rand('seed', 1);
  atx_bits = round(rand(1,code_param.data_bits_per_frame));
  tx_bits = []; tx_codewords = [];
  for f=1:interleave_frames
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
  Nerrs_coded_log = Nerrs_log = [];
  error_positions = [];
  Nerrs = Nerrs_coded = 0;
  
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

    % If looking for sync: check raw BER on frame just received
    % against all possible positions in the interleaver frame.

    % state machine(s) for modem and interleaver sync ------------------------------------

    if strcmp(states.sync_state,'searching') 
      [timing_valid states] = ofdm_sync_search(states, rxbuf_in);
    end
    
    if strcmp(states.sync_state,'synced') || strcmp(states.sync_state,'trial_sync')
      [rx_bits states aphase_est_pilot_log arx_np arx_amp] = ofdm_demod(states, rxbuf_in);
      rx_uw = rx_bits(1:states.Nuwbits);
      
      % we are in sync so log modem states

      rx_np_log = [rx_np_log arx_np];
      timing_est_log = [timing_est_log states.timing_est];
      delta_t_log = [delta_t_log states.delta_t];
      foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
      phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];

      % update sliding windows of rx-ed symbols and symbol amplitudes,
      % discarding UW and txt symbols at start of each modem frame

      rx_np(1:Nsymbolsperinterleavedframe-Nsymbolsperframe) = rx_np(Nsymbolsperframe+1:Nsymbolsperinterleavedframe);
      rx_np(Nsymbolsperinterleavedframe-Nsymbolsperframe+1:Nsymbolsperinterleavedframe) = arx_np(Nuwtxtsymbolsperframe+1:end);
      rx_amp(1:Nsymbolsperinterleavedframe-Nsymbolsperframe) = rx_amp(Nsymbolsperframe+1:Nsymbolsperinterleavedframe);
      rx_amp(Nsymbolsperinterleavedframe-Nsymbolsperframe+1:Nsymbolsperinterleavedframe) = arx_amp(Nuwtxtsymbolsperframe+1:end);

      % just single frame for now, so trival interlave sync
      % TODO: work out a way to get interleaver sync
      
      if (mod(frame_count,interleave_frames) == 0)
        
        % de-interleave QPSK symbols and symbol amplitudes

        arx_np = gp_deinterleave(rx_np);
        arx_amp = gp_deinterleave(rx_amp);

        % measure uncoded bit errors over interleaver frame

        rx_bits_raw = [];
        for s=1:Nsymbolsperinterleavedframe
          rx_bits_raw = [rx_bits_raw qpsk_demod(rx_np(s))];
        end
        for ff=1:interleave_frames
          st = (ff-1)*Ncodedbitsperframe+1; en = st+Ncodedbitsperframe-1;
          errors = xor(tx_bits_raw(st:en), rx_bits_raw(st:en));
          Nerrs = sum(errors);
          Terrs += Nerrs;
          Nerrs_log = [Nerrs_log Nerrs];
          Tbits += Ncodedbitsperframe;
        end
        
        % LDPC decode
        %  note: ldpc_errors can be used to measure raw BER
        %        std CML library doesn't have an indication of convergence

        rx_bits = [];
        for ff=1:interleave_frames
          st = (ff-1)*Ncodedbitsperframe/bps+1; en = st + Ncodedbitsperframe/bps - 1;
          [rx_codeword ldpc_errors] = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, arx_np(st:en), min(EsNo,30), arx_amp(st:en));
          rx_bits = [rx_bits rx_codeword(1:code_param.data_bits_per_frame)];          
        end

        errors_coded = xor(tx_bits, rx_bits);
        Nerrs_coded = sum(errors_coded);
        if Nerrs_coded < 0.2
          Terrs_coded += Nerrs_coded;
          Tbits_coded += code_param.data_bits_per_frame*interleave_frames;
          error_positions = [error_positions errors_coded];
        
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
            Nerrs_coded_log = [Nerrs_coded_log Nvocpacketerrs];
          end
        end
      end
 
    end
    
    states = sync_state_machine(states, rx_uw);

    if states.verbose
      printf("f: %2d state: %-10s uw_errors: %2d %1d Nerrs: %3d Nerrs_coded: %3d foff: %3.1f\n",
             f, states.last_sync_state, states.uw_errors, states.sync_counter, Nerrs, Nerrs_coded, states.foff_est_hz);
    end

    % act on any events returned by modem sync state machine
    
    if states.sync_start
      Nerrs_log = [];
      Terrs = Tbits = frame_count = 0;
      Tpacketerrs = Tpackets = 0;
      Terrs_coded = Tbits_coded = 0;
      error_positions = Nerrs_coded_log = [];
    end

  end

  printf("Coded BER: %5.4f Tbits: %5d Terrs: %5d PER: %5.4f Tpacketerrs: %5d Tpackets: %5d\n", 
         Terrs_coded/Tbits_coded, Tbits_coded, Terrs_coded, Tpacketerrs/Tpackets, Tpacketerrs, Tpackets);
  printf("Raw BER..: %5.4f Tbits: %5d Terrs: %5d\n", Terrs/Tbits, Tbits, Terrs);

  figure(1); clf; 
  plot(rx_np_log,'+');
  mx = max(abs(rx_np_log));
  axis([-mx mx -mx mx]);
  title('Scatter');

  figure(2); clf;
  plot(phase_est_pilot_log,'g+', 'markersize', 5); 
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
  mx = max(abs(foff_est_hz_log));
  axis([1 max(Nframes,2) -mx mx]);
  title('Fine Freq');
  ylabel('Hz')

  figure(5); clf;
  subplot(211)
  stem(Nerrs_log);
  title('Uncoded errrors/modem frame')
  axis([1 length(Nerrs_log) 0 Nbitsperframe*rate/2]);
  subplot(212)
  stem(Nerrs_coded_log);
  title('Coded errrors/vocoder frame')
  axis([1 length(Nerrs_coded_log) 0 Nbitspervocframe/2]);

  if nargin == 3
    fep = fopen(error_pattern_filename, "wb");
    fwrite(fep, error_positions, "short");
    fclose(fep);
  end
endfunction
