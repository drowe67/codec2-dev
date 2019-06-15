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

function time_to_sync = ofdm_ldpc_rx(filename, mode="700D", interleave_frames = 1, error_pattern_filename, start_secs, len_secs)
  ofdm_lib;
  ldpc;
  gp_interleaver;
  more off;

  % init modem

  [bps Rs Tcp Ns Nc] = ofdm_init_mode(mode);
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;
  states.verbose = 1;

  mod_order = 4; bps = 2; modulation = 'QPSK'; mapping = 'gray';
  demod_type = 0; decoder_type = 0; max_iterations = 100;

  EsNo = 3; % TODO: fixme
  printf("EsNo fixed at %f - need to est from channel\n", EsNo);
  
  % some constants used for assembling modem frames
  
  [code_param Nbitspercodecframe Ncodecframespermodemframe Nunprotect] = codec_to_frame_packing(states, mode);

  % load real samples from file

  Ascale= states.amp_scale/2.0;  % /2 as real signal has half amplitude
  frx=fopen(filename,"rb"); rx = fread(frx, Inf, "short")/Ascale; fclose(frx);
  length(rx)
  if (nargin >= 5) printf("start_secs: %d\n", start_secs); rx = rx(start_secs*Fs+1:end); end
  if (nargin >= 6) printf("len_secs: %d\n", len_secs); rx = rx(1:len_secs*Fs); end
  length(rx)
  Nsam = length(rx); Nframes = floor(Nsam/Nsamperframe);
  prx = 1;

  % OK generate tx frame for BER calcs
  %   We just use a single test frame of bits as it makes interleaver sync
  %   easier than using a test frame of bits that spans the entire interleaver
  %   frame.  Doesn't affect operation with the speech codec.
  
  if strcmp(mode, "700D")
    codec_bits = round(ofdm_rand(code_param.data_bits_per_frame)/32767);
  else
    codec_bits = round(ofdm_rand(Ncodecframespermodemframe*Nbitspercodecframe)/32767);
  end

  [frame_bits bits_per_frame] = assemble_frame(states, code_param, mode, codec_bits, Ncodecframespermodemframe, Nbitspercodecframe);

  % Some handy constants, "frame" refers to modem frame less UW and
  % txt bits.  For 700D this is comprised on one LDPC codeword, for
  % 2200 "codeword" is a mix of several LDPC codewords and some
  % unprotected bits.
  
  Ncodedbitsperframe = bits_per_frame;
  Nsymbolsperframe = bits_per_frame/bps;
  Nuwtxtsymbolsperframe = (Nuwbits+Ntxtbits)/bps;
  Nsymbolsperinterleavedframe = interleave_frames*Nsymbolsperframe;

  % buffers for interleaved frames

  rx_np = zeros(1, Nsymbolsperinterleavedframe);
  rx_amp = zeros(1, Nsymbolsperinterleavedframe);
  rx_uw = [];
  
  tx_bits = []; tx_frames = [];
  for f=1:interleave_frames
    tx_bits = [tx_bits codec_bits];
    tx_frames = [tx_frames frame_bits];
  end

  % used for rx frame sync on interleaved symbols - we demod the
  % entire interleaved frame to raw bits

  tx_symbols = [];
  for s=1:Nsymbolsperinterleavedframe
    tx_symbols = [tx_symbols qpsk_mod( tx_frames(2*(s-1)+1:2*s) )];
  end

  tx_symbols = gp_interleave(tx_symbols);

  tx_bits_raw = [];
  for s=1:Nsymbolsperinterleavedframe
    tx_bits_raw = [tx_bits_raw qpsk_demod(tx_symbols(s))];
  end

  % init logs and BER stats

  rx_bits = []; rx_np_log = []; timing_est_log = []; delta_t_log = []; foff_est_hz_log = [];
  phase_est_pilot_log = [];
  sig_var_log = []; noise_var_log = [];
  Terrs = Tbits = Terrs_coded = Tbits_coded = 0;
  Nerrs_coded_log = Nerrs_log = [];
  error_positions = [];
  Nerrs_coded = Nerrs_raw = zeros(1, interleave_frames);
  paritychecks = [0];
  time_to_sync = -1;
  
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

    % If looking for sync: check raw BER on frame just received
    % against all possible positions in the interleaver frame.

    % state machine(s) for modem and interleaver sync ------------------------------------

    if strcmp(states.sync_state,'search') 
      [timing_valid states] = ofdm_sync_search(states, rxbuf_in);
    end

    if strcmp(states.sync_state,'synced') || strcmp(states.sync_state,'trial')
      [rx_bits states aphase_est_pilot_log arx_np arx_amp] = ofdm_demod(states, rxbuf_in);
      [rx_uw payload_syms payload_amps txt_bits] = disassemble_modem_frame(states, arx_np, arx_amp);
           
      % we are in sync so log modem states

      rx_np_log = [rx_np_log arx_np];
      timing_est_log = [timing_est_log states.timing_est];
      delta_t_log = [delta_t_log states.delta_t];
      foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
      phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];
      sig_var_log = [sig_var_log states.sig_var];
      noise_var_log = [noise_var_log states.noise_var];
      
      % update sliding windows of rx-ed symbols and symbol amplitudes,
      % discarding UW and txt symbols at start of each modem frame

      rx_np(1:Nsymbolsperinterleavedframe-Nsymbolsperframe) = rx_np(Nsymbolsperframe+1:Nsymbolsperinterleavedframe);
      rx_np(Nsymbolsperinterleavedframe-Nsymbolsperframe+1:Nsymbolsperinterleavedframe) = payload_syms;
      rx_amp(1:Nsymbolsperinterleavedframe-Nsymbolsperframe) = rx_amp(Nsymbolsperframe+1:Nsymbolsperinterleavedframe);
      rx_amp(Nsymbolsperinterleavedframe-Nsymbolsperframe+1:Nsymbolsperinterleavedframe) = payload_amps;
           
      mean_amp = states.mean_amp;
      
      % de-interleave QPSK symbols and symbol amplitudes

      rx_np_de = gp_deinterleave(rx_np);
      rx_amp_de = gp_deinterleave(rx_amp);
      
      % Interleaver Sync:
      %   Needs to work on any data
      %   Use indication of LDPC convergence, may need to patch CML code for that
      %   Attempt a decode on every frame, when it converges we have sync

      next_sync_state_interleaver = states.sync_state_interleaver;

      if strcmp(states.sync_state_interleaver,'search')
        Nerrs = 0;
        if strcmp(mode, "700D")
          % using LDPC decoder to obtain interleaver sync only supported for 700D so far, lets assume 2200 only works
          % with interleave_frames = 1.
          st = 1; en = Ncodedbitsperframe/bps;
          [rx_codeword paritychecks] = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, rx_np_de(st:en)/mean_amp, min(EsNo,30), rx_amp_de(st:en)/mean_amp);
          Nerrs = code_param.data_bits_per_frame - max(paritychecks);
          %printf("Nerrs: %d\n", Nerrs);
        end
        
        % note we just go straight into sync if interleave_frames == 1
        
        if (Nerrs < 10) || (interleave_frames == 1)
          % sucessful(ish) decode!
          next_sync_state_interleaver = 'synced';
          states.frame_count_interleaver = interleave_frames;
        end
      end

      states.sync_state_interleaver = next_sync_state_interleaver;
            
      if strcmp(states.sync_state_interleaver,'synced') && (states.frame_count_interleaver == interleave_frames)
        states.frame_count_interleaver = 0;
        Nerrs_raw = Nerrs_coded = zeros(1, interleave_frames);

        %printf("decode!\n");
        
        % measure uncoded bit errors over interleaver frame

        rx_bits_raw = [];
        for s=1:Nsymbolsperinterleavedframe
          rx_bits_raw = [rx_bits_raw qpsk_demod(rx_np_de(s))];
        end
        for ff=1:interleave_frames
          st = (ff-1)*Ncodedbitsperframe+1; en = st+Ncodedbitsperframe-1;
          errors = xor(frame_bits, rx_bits_raw(st:en));
          Nerrs = sum(errors);
          Nerrs_log = [Nerrs_log Nerrs];
          Nerrs_raw(ff) += Nerrs;
          Tbits += Ncodedbitsperframe;
          Terrs += Nerrs;
        end
        
        % LDPC decode
        %  note: ldpc_errors can be used to measure raw BER
        %        std CML library doesn't have an indication of convergence

        rx_bits = [];
        for ff=1:interleave_frames

          % this is where we need to modify things for "2200".  LDPC decode then rebuild codec_bits.  Hmm, would be easier
          % with a custom designed LDPC code.

          if strcmp(mode, "700D")
            st = (ff-1)*Nsymbolsperframe+1; en = st + Nsymbolsperframe-1;
            [rx_codeword paritychecks] = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, rx_np_de(st:en)/mean_amp, min(EsNo,30), rx_amp_de(st:en)/mean_amp);
            arx_bits = rx_codeword(1:code_param.data_bits_per_frame);
            errors = xor(codec_bits, arx_bits);
            Nerrs  = sum(errors);
            Tbits_coded += code_param.data_bits_per_frame;
          else
            % "2200" mode, run LDPC decoder twice, re-assemble codec bit stream
            
            assert(interleave_frames == 1);
            Ncodewordsymbolsperframe = code_param.code_bits_per_frame/bps;
            data_bits_per_frame = code_param.data_bits_per_frame;
            code_bits_per_frame = code_param.code_bits_per_frame;
            st = 1; en = st + Ncodewordsymbolsperframe - 1;
            rx_codeword1 = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, rx_np_de(st:en)/mean_amp, min(EsNo,30), rx_amp_de(st:en)/mean_amp);
            st += Ncodewordsymbolsperframe; en += Ncodewordsymbolsperframe;
            rx_codeword2 = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, rx_np_de(st:en)/mean_amp, min(EsNo,30), rx_amp_de(st:en)/mean_amp);

            protected_bits = [rx_codeword1(1:data_bits_per_frame) rx_codeword2(1:data_bits_per_frame)];
            unprotected_bits = rx_bits_raw(2*code_bits_per_frame+1:2*code_bits_per_frame+Nunprotect);
            %printf("lprot: %d lunprot: %d\n", length(protected_bits), length(unprotected_bits));
            
            % OK now reconstruct codec frames from protected and unprotected bits

            Nprotectedbitspercodecframe = 2*code_param.data_bits_per_frame/Ncodecframespermodemframe;
            Nunprotectedbitspercodecframe = Nbitspercodecframe - Nprotectedbitspercodecframe;
            rx_codec_bits = []; Nerrs = 0;
            for i=1:Ncodecframespermodemframe
              a = (i-1)*Nprotectedbitspercodecframe + 1;
              b = a + Nprotectedbitspercodecframe-1;
              c = (i-1)*Nunprotectedbitspercodecframe+1;
              d = c + Nunprotectedbitspercodecframe-1;
              rx_codec_bits = [rx_codec_bits protected_bits(a:b) unprotected_bits(c:d)];

              % note we are just counting errors in coded part
              
              e  = (i-1)*Nbitspercodecframe + 1;
              g  = e + Nprotectedbitspercodecframe-1;
              Nerrs += sum(xor(codec_bits(e:g), protected_bits(a:b)));
              Tbits_coded += Nprotectedbitspercodecframe;
            end
            arx_bits = rx_codec_bits;

            % to generate codec error bit stream compare reconstructed rx codec frames to original
            
            codec_errors = xor(codec_bits, rx_codec_bits);
            error_positions = [error_positions codec_errors]; 
          end

          rx_bits = [rx_bits arx_bits];
          
          Nerrs_coded(ff) = Nerrs;
          Terrs_coded += Nerrs;
          Nerrs_coded_log = [Nerrs_coded_log Nerrs];
        end
      end
    end
    
    states = sync_state_machine(states, rx_uw);

    if states.verbose
      r = mod(states.frame_count_interleaver,  interleave_frames)+1;
      pcc = max(paritychecks);
      iter = 0;
      for i=1:length(paritychecks)
        if paritychecks(i) iter=i; end
      end
      printf("f: %3d st: %-6s euw: %2d %1d ist: %-6s eraw: %3d ecdd: %3d iter: %3d pcc: %3d foff: %4.1f\n",
             f, states.last_sync_state, states.uw_errors, states.sync_counter, states.last_sync_state_interleaver,
             Nerrs_raw(r), Nerrs_coded(r), iter, pcc, states.foff_est_hz);
      % detect a sucessful sync
      if (time_to_sync < 0) && (strcmp(states.sync_state,'synced') || strcmp(states.sync_state,'trial'))
        if (pcc > 80) && (iter != 100)
          time_to_sync = f*Nsamperframe/Fs;
        end
      end
    end

    % act on any events returned by modem sync state machine
    
    if states.sync_start
      Nerrs_raw = Nerrs_coded = zeros(1, interleave_frames);
      Nerrs_log = [];
      Terrs = Tbits = 0;
      Tpacketerrs = Tpackets = 0;
      Terrs_coded = Tbits_coded = 0;
      error_positions = Nerrs_coded_log = [];
    end

  end

  printf("Raw BER..: %5.4f Tbits: %5d Terrs: %5d\n", Terrs/(Tbits+1E-12), Tbits, Terrs);
  printf("Coded BER: %5.4f Tbits: %5d Terrs: %5d\n", Terrs_coded/(Tbits_coded+1E-12), Tbits_coded, Terrs_coded);

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

  if length(Nerrs_log)
    figure(5); clf;
    subplot(211)
    stem(Nerrs_log);
    title('Uncoded errrors/modem frame')
    axis([1 length(Nerrs_log) 0 Nbitsperframe*rate/2]);
    if length(Nerrs_coded_log)
      subplot(212)
      stem(Nerrs_coded_log);
      title('Coded errors/mode frame')
      axis([1 length(Nerrs_coded_log) 0 Nbitsperframe/2]);
    end
  end
  
  figure(6); clf;
  snr_estdB = 10*log10(sig_var_log) - 10*log10(noise_var_log) + 10*log10(Nc*Rs/3000);
  snr_smoothed_estdB = filter(0.1,[1 -0.9],snr_estdB);
  plot(snr_smoothed_estdB);
  title('Signal and Noise Power estimates');
  ylabel('SNR (dB)')
  
  if (nargin == 4) && strlen(error_pattern_filename)
    fep = fopen(error_pattern_filename, "wb");
    fwrite(fep, error_positions, "uchar");
    fclose(fep);
  end
endfunction
