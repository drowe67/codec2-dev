% ofdm_rx.m
% David Rowe May 2018
%
% OFDM file based rx to unit test core OFDM modem.  See also
% ofdm_ldpc_rx which includes LDPC and interleaving, and ofdm_demod.c

#{
  TODO:
    [X] single frame based sync state machine
        + that doesn't depend on payload data
    [ ] make robust to fading
        + what are win conditions?
        + when to resync?
        + how long to hang on?
    [ ] bring into line with C version, e.g. how test bits are fed in
    [ ] audio levels in 16 bit shorts compatable with other modes
#}

function ofdm_rx(filename, error_pattern_filename)
  ofdm_lib;
  more off;

  % init modem

  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 16; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;
  states.verbose = 1;

  % load real samples from file

  Ascale= 2E5*1.1491;
  frx=fopen(filename,"rb"); rx = 2*fread(frx, Inf, "short")/4E5; fclose(frx);
  Nsam = length(rx); Nframes = floor(Nsam/Nsamperframe);
  prx = 1;

  % OK re-generate tx frame for BER calcs

  rand('seed', 1);
  tx_bits = round(rand(1,Nbitsperframe));

  % init logs and BER stats

  rx_bits = []; rx_np_log = []; timing_est_log = []; delta_t_log = []; foff_est_hz_log = [];
  phase_est_pilot_log = [];
  Terrs = Tbits = Terrs_coded = Tbits_coded = Tpackets = Tpacketerrs = 0;
  Nbitspervocframe = 28;
  Nerrs_coded_log = Nerrs_log = [];
  error_positions = [];

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

    [rx_bits states aphase_est_pilot_log arx_np arx_amp] = ofdm_demod(states, rxbuf_in);

    errors = xor(tx_bits, rx_bits);
    Nerrs = sum(errors);
    aber = Nerrs/Nbitsperframe;
    
    frame_count++;

    printf("f: %d state: %s Nerrs: %d aber: %3.2f\n", f, state, Nerrs, aber);

    % If looking for sync: check raw BER on frame just received
    % against all possible positions in the interleaver frame.

    % iterate state machine ------------------------------------

    next_state = state;
    if strcmp(state,'searching')  

      % If looking for sync: check raw BER on frame just received
      % against all possible positions in the interleaver frame.

      if aber < 0.1
        next_state = 'synced';
        % make sure we get an interleave frame with correct freq offset
        % note this introduces a lot of delay, a better idea would be to
        % run demod again from interleave_frames back with now-known freq offset
      end
    end

    if strcmp(state,'synced')  
      if Nerrs/Nbitsperframe > 0.2
        %next_state = 'searching';
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
    
    if strcmp(state,'synced')

      % we are in sync so log states

      rx_np_log = [rx_np_log arx_np];
      timing_est_log = [timing_est_log states.timing_est];
      delta_t_log = [delta_t_log states.delta_t];
      foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
      phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];

      % measure uncoded bit errors on modem frame

      if aber < 0.2
        Terrs += Nerrs;
        Nerrs_log = [Nerrs_log Nerrs];
        Tbits += Nbitsperframe;
      end
    end
  end

  printf("BER..: %5.4f Tbits: %5d Terrs: %5d\n", Terrs/Tbits, Tbits, Terrs);

  figure(1); clf; 
  plot(rx_np_log,'+');
  mx = max(abs(rx_np_log));
  %axis([-mx mx -mx mx]);
  title('Scatter');

  figure(2); clf;
  plot(phase_est_pilot_log(:,2:Nc),'g+', 'markersize', 5); 
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
  mx = max(abs(foff_est_hz_log))+1;
  axis([1 max(Nframes,2) -mx mx]);
  title('Fine Freq');
  ylabel('Hz')

  figure(5); clf;
  stem(Nerrs_log);
  title('Errors/modem frame')
  axis([1 length(Nerrs_log) 0 Nbitsperframe*rate/2]);

  if nargin == 3
    fep = fopen(error_pattern_filename, "wb");
    fwrite(fep, error_positions, "short");
    fclose(fep);
  end
endfunction
