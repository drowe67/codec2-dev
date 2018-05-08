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

  EsNo = 3; % TODO: fixme
  printf("EsNo fixed at %f - need to est from channel\n", EsNo);
  
  init_cml('~/cml/');
  load HRA_112_112.txt
  [code_param framesize rate] = ldpc_init_user(HRA_112_112, modulation, mod_order, mapping);
  assert(Nbitsperframe == (code_param.code_bits_per_frame + states.Nuwbits + states.Ntxtbits));

  % load real samples from file

  Ascale= states.amp_scale/2.0;  % /2 as real signal has half amplitude
  frx=fopen(filename,"rb"); rx = fread(frx, Inf, "short")/Ascale; fclose(frx);
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
  %   We just use a single test frame of bits as it makes interleaver sync
  %   easier than using a test frame of bits that spans the entire interleaver
  %   frame.  Doesn't affect operation with the speech codec operation.
  
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
  sig_var_log = []; noise_var_log = [];
  Terrs = Tbits = Terrs_coded = Tbits_coded = Tpackets = Tpacketerrs = 0;
  Nbitspervocframe = 28;
  Nerrs_coded_log = Nerrs_log = [];
  error_positions = [];
  Nerrs_coded = Nerrs_raw = zeros(1, interleave_frames);

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
        st = 1; en = Ncodedbitsperframe/bps;
        [rx_codeword parity_checks] = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, rx_np_de(st:en)/mean_amp, min(EsNo,30), rx_amp_de(st:en)/mean_amp);
        Nerrs = code_param.data_bits_per_frame - max(parity_checks);
        %printf("Nerrs: %d\n", Nerrs);

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
          errors = xor(acodeword, rx_bits_raw(st:en));
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
          st = (ff-1)*Ncodedbitsperframe/bps+1; en = st + Ncodedbitsperframe/bps - 1;
          [rx_codeword ldpc_errors] = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, rx_np_de(st:en)/mean_amp, min(EsNo,30), rx_amp_de(st:en)/mean_amp);
          rx_bits = [rx_bits rx_codeword(1:code_param.data_bits_per_frame)];
          errors = xor(atx_bits, rx_codeword(1:code_param.data_bits_per_frame));
          Nerrs  = sum(errors);
          Nerrs_coded(ff) = Nerrs;
          Terrs_coded += Nerrs;
          Tbits_coded += code_param.data_bits_per_frame;       
        end

        % additional measure - PER on vocoder frames

        errors_coded = xor(tx_bits, rx_bits);
        Nerrs = sum(errors_coded);
        if Nerrs/(code_param.data_bits_per_frame*interleave_frames) < 0.2
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
      r = mod(states.frame_count_interleaver,  interleave_frames)+1;
      printf("f: %3d st: %-6s uw_errs: %2d %1d inter_st: %-6s inter_fr: %2d Nerrs_raw: %3d Nerrs_coded: %3d foff: %4.1f\n",
             f, states.last_sync_state, states.uw_errors, states.sync_counter, states.last_sync_state_interleaver, states.frame_count_interleaver,
             Nerrs_raw(r), Nerrs_coded(r), states.foff_est_hz);
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
  printf("Codec PER: %5.4f Tpkts: %5d Terrs: %5d\n", Tpacketerrs/(Tpackets+1E-12), Tpackets, Tpacketerrs);

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
      title('Coded errrors/vocoder frame')
      axis([1 length(Nerrs_coded_log) 0 Nbitspervocframe/2]);
    end
  end
  
  figure(6); clf;
  snr_estdB = 10*log10(sig_var_log) - 10*log10(noise_var_log) + 10*log10(Nc*Rs/3000);
  snr_smoothed_estdB = filter(0.1,[1 -0.9],snr_estdB);
  plot(snr_smoothed_estdB);
  title('Signal and Noise Power estimates');
  ylabel('SNR (dB)')
  
  if nargin == 3
    fep = fopen(error_pattern_filename, "wb");
    fwrite(fep, error_positions, "short");
    fclose(fep);
  end
endfunction
