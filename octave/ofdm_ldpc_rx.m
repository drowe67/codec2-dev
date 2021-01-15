% ofdm_ldpc_rx.m
% David Rowe April 2017
%
% OFDM file based rx, with LDPC and interleaver, Octave version of src/ofdm_demod.c

function time_to_sync = ofdm_ldpc_rx(filename, mode="700D", error_pattern_filename, start_secs, len_secs=0)
  ofdm_lib;
  ldpc;
  gp_interleaver;
  more off;
  pkg load signal;

  % optional DPSK testing
  dpsk = 0;
  if strcmp(mode,"700D-DPSK")
    mode = "700D"; dpsk = 1;
  end
  if strcmp(mode,"2020-DPSK")
    mode = "2020"; dpsk = 1;
  end

  % init modem

  config = ofdm_init_mode(mode);
  states = ofdm_init(config);
  ofdm_load_const;
  states.verbose = 1;
  states.dpsk = dpsk;

  % some constants used for assembling modem frames

  [code_param Nbitspercodecframe Ncodecframespermodemframe] = codec_to_frame_packing(states, mode);

  % load real samples from file

  Ascale= states.amp_scale/2.0;  % /2 as real signal has half amplitude
  frx=fopen(filename,"rb"); rx = fread(frx, Inf, "short")/Ascale; fclose(frx);
  if (nargin >= 5) printf("start_secs: %d\n", start_secs); rx = rx(start_secs*Fs+1:end); end
  if (nargin >= 6) printf("len_secs: %d\n", len_secs); rx = rx(1:len_secs*Fs); end
  Nsam = length(rx); Nframes = floor(Nsam/Nsamperframe);
  prx = 1;

  % Generate tx frame for BER calcs

  payload_bits = round(ofdm_rand(code_param.data_bits_per_frame)/32767);
  tx_bits = fec_encode(states, code_param, mode, payload_bits, Ncodecframespermodemframe, Nbitspercodecframe);

  % Some handy constants

  Nsymsperframe = Nbitsperframe/bps;
  Nsymsperpacket = Nbitsperpacket/bps;
  Ncodedbitsperpacket = code_param.coded_bits_per_frame;
  Ncodedsymsperpacket = code_param.coded_syms_per_frame;

  % init logs and BER stats

  rx_bits = []; rx_np_log = []; timing_est_log = []; delta_t_log = []; foff_est_hz_log = [];
  channel_est_pilot_log = [];
  sig_var_log = []; noise_var_log = []; EsNo_log = []; mean_amp_log = [];
  Terrs = Tbits = Terrs_coded = Tbits_coded = Perrs_coded = 0;
  Nerrs_coded_log = Nerrs_log = [];
  error_positions = [];
  Nerrs_coded = Nerrs_raw = 0;
  paritychecks = [0];
  time_to_sync = -1; EsNo = 1;
  rx_uw = zeros(1,states.Nuwbits);

  rx_syms = zeros(1,Nsymsperpacket); rx_amps = zeros(1,Nsymsperpacket);
  packet_count = frame_count = 0;

  #{
  % 'prime' rx buf to get correct coarse timing (for now)

  prx = 1;
  nin = Nsamperframe+2*(M+Ncp);
  states.rxbuf(Nrxbuf-nin+1:Nrxbuf) = rx(prx:nin);
  prx += nin;
  #}

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

    if states.verbose
      printf("f: %3d nin: %4d st: %-6s ", f, states.nin, states.sync_state);
    end

    if strcmp(states.sync_state,'search')
      [timing_valid states] = ofdm_sync_search(states, rxbuf_in);
    end

    if strcmp(states.sync_state,'synced') || strcmp(states.sync_state,'trial')
      % accumulate a buffer of data symbols for this packet
      rx_syms(1:end-Nsymsperframe) = rx_syms(Nsymsperframe+1:end);
      rx_amps(1:end-Nsymsperframe) = rx_amps(Nsymsperframe+1:end);
      [states rx_bits achannel_est_pilot_log arx_np arx_amp] = ofdm_demod(states, rxbuf_in);
      rx_syms(end-Nsymsperframe+1:end) = arx_np;
      rx_amps(end-Nsymsperframe+1:end) = arx_amp;

      rx_uw = extract_uw(states, rx_syms(end-Nuwframes*Nsymsperframe+1:end), rx_amps(end-Nuwframes*Nsymsperframe+1:end));

      % We need the full packet of symbols before disassembling and checking for bit errors
      if (states.modem_frame == (states.Np-1))
        packet_count++;

        % unpack, de-interleave PSK symbols and symbol amplitudes
        [rx_uw_unused payload_syms payload_amps txt_bits] = disassemble_modem_packet(states, rx_syms, rx_amps);
        payload_syms_de = gp_deinterleave(payload_syms);
        payload_amps_de = gp_deinterleave(payload_amps);

        % Count uncoded (raw) errors
        rx_bits = zeros(1,Ncodedbitsperpacket);
        for s=1:Ncodedsymsperpacket
          if bps == 2 rx_bits(2*s-1:2*s) = qpsk_demod(payload_syms_de(s)); end
          if bps == 4 rx_bits(bps*(s-1)+1:bps*s) = qam16_demod(states.qam16,payload_syms_de(s), payload_amps_de(s)); end
        end
        errors = xor(tx_bits, rx_bits);
        Nerrs = sum(errors);
        Nerrs_log = [Nerrs_log Nerrs]; Nerrs_raw = Nerrs;
        Terrs += Nerrs;
        Tbits += Nbitsperpacket;

        % LDPC decode

        % keep earlier mean amplitude estimator for compatability with 700D
        if states.amp_est_mode == 0
          mean_amp = states.mean_amp;
        else
          mean_amp = mean(payload_amps_de)+1E-12;
        end
        mean_amp_log = [mean_amp_log mean_amp];

        % used fixed EsNo est, as EsNo estimator for QAM not working very well at this stage
        EsNo = 10^(states.EsNodB/10);

        % TODO 2020 support for padding with known data bits

        [rx_codeword paritychecks] = ldpc_dec(code_param, mx_iter=100, demod=0, dec=0, ...
                                              payload_syms_de/mean_amp, EsNo, payload_amps_de/mean_amp);
        rx_bits = rx_codeword(1:code_param.data_bits_per_frame);
        errors = xor(payload_bits, rx_bits);
        Nerrs_coded  = sum(errors);
        EsNo_log = [EsNo_log EsNo];

        % reset EsNo estimator for next packet
        states.sum_sig_var = states.sum_noise_var = 0;

        if Nerrs_coded Perrs_coded++; end
        Terrs_coded += Nerrs_coded;
        Tbits_coded += code_param.data_bits_per_frame;
        Nerrs_coded_log = [Nerrs_coded_log Nerrs_coded];
      end

      % we are in sync so log modem states

      rx_np_log = [rx_np_log arx_np];
      timing_est_log = [timing_est_log states.timing_est];
      delta_t_log = [delta_t_log states.delta_t];
      foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
      channel_est_pilot_log = [channel_est_pilot_log; achannel_est_pilot_log];
      sig_var_log = [sig_var_log states.sig_var];
      noise_var_log = [noise_var_log states.noise_var];
      frame_count++;
    end

    states = sync_state_machine(states, rx_uw);

    if states.verbose
      if strcmp(states.last_sync_state,'synced') || strcmp(states.last_sync_state,'trial')
        pcc = max(paritychecks);
        iter = 0;
        for i=1:length(paritychecks)
          if paritychecks(i) iter=i; end
        end
        % complete logging line
        printf("euw: %3d %d mf: %2d pbw: %s eraw: %3d ecod: %3d iter: %3d pcc: %3d EsNo: %4.2f foff: %4.1f",
               states.uw_errors, states.sync_counter, states.modem_frame, states.phase_est_bandwidth(1),
               Nerrs_raw, Nerrs_coded, iter, pcc, 10*log10(EsNo), states.foff_est_hz);
        % detect a sucessful sync (for tests calling this function)
        if (time_to_sync < 0) && (strcmp(states.sync_state,'synced') || strcmp(states.sync_state,'trial'))
          if (pcc > 80) && (iter != 100)
            time_to_sync = f*Nsamperframe/Fs;
          end
        end
      end
      printf("\n");
    end

    % act on any events returned by modem sync state machine

    if states.sync_start
      Nerrs_raw = Nerrs_coded = 0;
      Nerrs_log = [];
      Terrs = Tbits = 0;
      Tpacketerrs = Tpackets = 0;
      Terrs_coded = Tbits_coded = 0;
      error_positions = Nerrs_coded_log = [];
    end
  end

  printf("Raw BER..: %5.4f Tbits: %5d Terrs: %5d\n", Terrs/(Tbits+1E-12), Tbits, Terrs);
  printf("Coded BER: %5.4f Tbits: %5d Terrs: %5d\n", Terrs_coded/(Tbits_coded+1E-12), Tbits_coded, Terrs_coded);
  printf("Coded PER: %5.4f Pckts: %5d Perrs: %5d\n", Perrs_coded/(packet_count+1E-12), packet_count, Perrs_coded);

  if length(rx_np_log)
      figure(1); clf;
      plot(exp(j*pi/4)*rx_np_log(floor(end/4):floor(end-end/8)),'+');
      mx = 2*mean(abs(channel_est_pilot_log(:)));
      axis([-mx mx -mx mx]);
      title('Scatter');

      figure(2); clf;
      plot(angle(channel_est_pilot_log),'g+', 'markersize', 5);
      title('Phase est');
      axis([1 length(channel_est_pilot_log) -pi pi]);

      figure(3); clf;
      amp_est = abs(channel_est_pilot_log);
      plot(amp_est,'g+', 'markersize', 5);
      title('Amp est');
      axis([1 length(channel_est_pilot_log) min(amp_est(:)) max(amp_est(:))]);

      figure(4); clf;
      subplot(211); plot(10*log10(EsNo_log+1E-12)); ylabel('EsNodB');
      subplot(212); plot(mean_amp_log); ylabel('mean amp');

      figure(5); clf;
      subplot(211)
      stem(delta_t_log)
      title('delta t');
      subplot(212)
      plot(timing_est_log);
      title('timing est');

      figure(6); clf;
      plot(foff_est_hz_log)
      mx = max(max(abs(foff_est_hz_log)),1);
      axis([1 max(Nframes,2) -mx mx]);
      title('Fine Freq');
      ylabel('Hz')
  end

  if length(Nerrs_log) > 1
    figure(7); clf;
    subplot(211)
    stem(Nerrs_log);
    title('Uncoded errrors/modem frame')
    axis([1 length(Nerrs_log) 0 Nbitsperpacket*0.2]);
    if length(Nerrs_coded_log)
      subplot(212)
      stem(Nerrs_coded_log);
      title('Coded errors/mode frame')
      axis([1 length(Nerrs_coded_log) 0 Nbitsperpacket*0.2]);
    end
  end

  figure(8); clf;
  snr_estdB = 10*log10(sig_var_log) - 10*log10(noise_var_log+1E-12) + 10*log10(Nc*Rs/3000);
  snr_smoothed_estdB = filter(0.1,[1 -0.9],snr_estdB);
  plot(snr_smoothed_estdB);
  title('Signal and Noise Power estimates');
  ylabel('SNR (dB)')

  figure(9); clf; plot_specgram(rx);
  if len_secs
    axis([0 len_secs 500 2500])
  end

  if (nargin == 4) && strlen(error_pattern_filename)
    fep = fopen(error_pattern_filename, "wb");
    fwrite(fep, error_positions, "uchar");
    fclose(fep);
  end
endfunction
