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

function ofdm_rx(filename)
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

  % fixed test frame of tx bits

  rand('seed', 100);
  tx_bits = round(rand(1,code_param.data_bits_per_frame));
  tx_codeword = LdpcEncode(tx_bits, code_param.H_rows, code_param.P_matrix);

  % init logs and BER stats

  rx_bits = []; rx_np = []; timing_est_log = []; delta_t_log = []; foff_est_hz_log = [];
  phase_est_pilot_log = [];
  Terrs = Tbits = 0;

  % load real samples from file

  Ascale= 2E5;
  frx=fopen(filename,"rb"); rx = 2*fread(frx, Inf, "short")/4E5; fclose(frx);
  Nsam = length(rx); Nframes = floor(Nsam/Nsamperframe);
  prx = 1;

  % 'prime' rx buf to get correct coarse timing (for now)

  prx = 1;
  nin = Nsamperframe+2*(M+Ncp);
  states.rxbuf(Nrxbuf-nin+1:Nrxbuf) = rx(prx:nin);
  prx += nin;

  state = 'searching';

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
    rx_codeword = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, arx_np, min(EsNo,30), arx_amp);
    rx_bits = rx_codeword(1:code_param.data_bits_per_frame);

    % measure errors and iterate start machine

    errors = xor(tx_bits, rx_bits);
    Nerrs = sum(errors);
    next_state = state;

    if strcmp(state,'searching')  
      if Nerrs == 0
        next_state = 'synced';
      end
    end
    state = next_state;

    if strcmp(state,'searching') 

      % attempt coarse timing estimate (i.e. detect start of frame)

      st = M+Ncp + Nsamperframe + 1; en = st + 2*Nsamperframe; 
      [ct_est foff_est] = coarse_sync(states, states.rxbuf(st:en), states.rate_fs_pilot_samples);
      if states.verbose
        printf("Nerrs: %d ct_est: %4d foff_est: %3.1f\n", Nerrs, ct_est, foff_est);
      end

      % calculate number of samples we need on next buffer to get into sync
     
      states.nin = Nsamperframe + ct_est - 1;

      % reset modem states

      states.sample_point = states.timing_est = 1;
      states.foff_est_hz = foff_est;

    else

      % we are in sync so log states and bit errors

      rx_np = [rx_np arx_np];
      timing_est_log = [timing_est_log states.timing_est];
      delta_t_log = [delta_t_log states.delta_t];
      foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
      phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];

      % measure bit errors

      Terrs += Nerrs;
      Nerrs_log(f) = Nerrs;
      Tbits += code_param.data_bits_per_frame;
    end
  end

  printf("Coded BER: %5.4f Tbits: %d Terrs: %d\n", Terrs/Tbits, Tbits, Terrs);

  figure(1); clf; 
  plot(rx_np,'+');
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
  mx = max(abs(foff_est_hz_log));
  axis([1 max(Nframes,2) -mx mx]);
  title('Fine Freq');

  figure(5); clf;
  %plot(Nerrs_log);
  
endfunction
