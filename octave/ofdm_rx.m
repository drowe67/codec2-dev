% ofdm_rx.m
% David Rowe May 2018
%
% OFDM file based rx to unit test core OFDM modem.  See also
% ofdm_ldpc_rx which includes LDPC and interleaving, and ofdm_demod.c


function ofdm_rx(filename, error_pattern_filename)
  ofdm_lib;
  more off;

  % init modem

  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 17; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;
  states.verbose = 0;

  % load real samples from file

  Ascale= 2E5*1.1491/2;
  frx=fopen(filename,"rb"); rx = fread(frx, Inf, "short")/Ascale; fclose(frx);
  Nsam = length(rx); Nframes = floor(Nsam/Nsamperframe);
  prx = 1;

  % OK re-generate tx frame for BER calcs

  rand('seed', 1);
  tx_bits = round(rand(1,Nbitsperframe));
  tx_bits(1:states.uw_len) = 0;   % insert UW
 
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
  %states.rxbuf(Nrxbuf-nin+1:Nrxbuf) = rx(prx:nin);
  %prx += nin;

  states.sync_state = states.last_sync_state = 'searching';
  states.uw_errors = 0;
  states.sync_counter = 0;
  states.sync_frame_count = 0;
  states.sync_start = 0;
  states.sync_end = 0;
  
  states.verbose = 1;

  Nerrs = 0; rx_uw = zeros(1,states.uw_len);
  
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
 
    if strcmp(states.sync_state,'searching') 
      [timing_valid states] = ofdm_sync_search(states, rxbuf_in);
    end
    
    if strcmp(states.sync_state,'synced') || strcmp(states.sync_state,'trial_sync')
      [rx_bits states aphase_est_pilot_log arx_np arx_amp] = ofdm_demod(states, rxbuf_in);
      rx_uw = rx_bits(1:states.uw_len);
      
      errors = xor(tx_bits, rx_bits);
      Nerrs = sum(errors);
      aber = Nerrs/Nbitsperframe;
    
      % we are in sync so log states

      rx_np_log = [rx_np_log arx_np];
      timing_est_log = [timing_est_log states.timing_est];
      delta_t_log = [delta_t_log states.delta_t];
      foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
      phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];

      % measure uncoded bit errors on modem frame

      Nerrs_log = [Nerrs_log Nerrs];
      Terrs += Nerrs;
      Tbits += Nbitsperframe;

      frame_count++;
    end
    
    states = sync_state_machine(states, rx_uw);

    if states.verbose
      printf("f: %2d state: %-10s uw_errors: %2d %1d Nerrs: %3d foff: %3.1f\n",
             f, states.last_sync_state, states.uw_errors, states.sync_counter, Nerrs, states.foff_est_hz);
    end

    % act on any events returned by state machine
    
    if states.sync_start
      Nerrs_log = [];
      Terrs = Tbits = frame_count = 0;
    end
  end

  printf("\nBER..: %5.4f Tbits: %5d Terrs: %5d\n", Terrs/Tbits, Tbits, Terrs);

  % If we have enough frames, calc BER discarding first few frames where freq
  % offset is adjusting

  Ndiscard = 20;
  if frame_count > Ndiscard
    Terrs -= sum(Nerrs_log(1:Ndiscard)); Tbits -= Ndiscard;
    printf("BER2.: %5.4f Tbits: %5d Terrs: %5d\n", Terrs/Tbits, Tbits, Terrs);
  end
  
  figure(1); clf; 
  plot(rx_np_log,'+');
  mx = 2*max(abs(rx_np_log));
  axis([-mx mx -mx mx]);
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
