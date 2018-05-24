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

  Ascale = states.amp_scale/2; % as input is a real valued signal
  frx=fopen(filename,"rb"); rx = fread(frx, Inf, "short")/Ascale; fclose(frx);
  Nsam = length(rx); Nframes = floor(Nsam/Nsamperframe);
  prx = 1;

  % OK re-generate tx frame for BER calcs

  tx_bits = create_ldpc_test_frame;

  % init logs and BER stats

  rx_bits = []; rx_np_log = []; timing_est_log = []; delta_t_log = []; foff_est_hz_log = [];
  phase_est_pilot_log = []; sig_var_log = []; noise_var_log = [];
  Terrs = Tbits = Terrs_coded = Tbits_coded = Tpackets = Tpacketerrs = frame_count = 0;
  Nbitspervocframe = 28;
  Nerrs_coded_log = Nerrs_log = [];
  error_positions = [];

  % 'prime' rx buf to get correct coarse timing (for now)

  prx = 1;
  nin = Nsamperframe+2*(M+Ncp);
  %states.rxbuf(Nrxbuf-nin+1:Nrxbuf) = rx(prx:nin);
  %prx += nin;
  
  states.verbose = 1;

  Nerrs = 0; rx_uw = zeros(1,states.Nuwbits);
  
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
 
    if strcmp(states.sync_state,'search') 
      [timing_valid states] = ofdm_sync_search(states, rxbuf_in);
    end
    
    if strcmp(states.sync_state,'synced') || strcmp(states.sync_state,'trial')
      [rx_bits states aphase_est_pilot_log arx_np arx_amp] = ofdm_demod(states, rxbuf_in);
      [rx_uw payload_syms payload_amps txt_bits] = disassemble_modem_frame(states, arx_np, arx_amp);
      
      errors = xor(tx_bits, rx_bits);
      Nerrs = sum(errors);
      aber = Nerrs/Nbitsperframe;
    
      % we are in sync so log states

      rx_np_log = [rx_np_log arx_np];
      timing_est_log = [timing_est_log states.timing_est];
      delta_t_log = [delta_t_log states.delta_t];
      foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
      phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];
      sig_var_log = [sig_var_log states.sig_var];
      noise_var_log = [noise_var_log states.noise_var];

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
      rx_np_log = [];
      sig_var_log = []; noise_var_log = [];
    end
  end

  printf("\nBER..: %5.4f Tbits: %5d Terrs: %5d\n", Terrs/(Tbits+1E-12), Tbits, Terrs);

  % If we have enough frames, calc BER discarding first few frames where freq
  % offset is adjusting

  Ndiscard = 20;
  if frame_count > Ndiscard
    Terrs -= sum(Nerrs_log(1:Ndiscard)); Tbits -= Ndiscard*Nbitsperframe;
    printf("BER2.: %5.4f Tbits: %5d Terrs: %5d\n", Terrs/Tbits, Tbits, Terrs);
  end

  %EsNo_est = mean(sig_var_log(floor(end/2):end))/mean(noise_var_log(floor(end/2):end));
  EsNo_est = mean(sig_var_log)/mean(noise_var_log);
  EsNo_estdB = 10*log10(EsNo_est);
  SNR_estdB = EsNo_estdB + 10*log10(Nc*Rs/3000);
  printf("Es/No est dB: % -4.1f SNR3k: %3.2f %f %f\n", EsNo_estdB, SNR_estdB, mean(sig_var_log), mean(noise_var_log));
  
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

  figure(6); clf;
  plot(10*log10(sig_var_log),'b;Es;');
  hold on;
  plot(10*log10(noise_var_log),'r;No;');
  snr_estdB = 10*log10(sig_var_log) - 10*log10(noise_var_log) + 10*log10(Nc*Rs/3000);
  snr_smoothed_estdB = filter(0.1,[1 -0.9],snr_estdB);
  plot(snr_smoothed_estdB,'g;SNR3k;');
  hold off;
  title('Signal and Noise Power estimates');

  if nargin == 2
    fep = fopen(error_pattern_filename, "wb");
    fwrite(fep, error_positions, "short");
    fclose(fep);
  end
endfunction
