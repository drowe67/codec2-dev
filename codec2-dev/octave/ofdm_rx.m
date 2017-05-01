% ofdm_rx.m
% David Rowe April 2017
%
% OFDM file based rx.


function ofdm_rx(filename)
  ofdm_lib;
  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 16; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;

  % fixed test frame of tx bits

  rand('seed', 100);
  tx_bits = rand(1,Nbitsperframe) > 0.5;

  % init logs and BER stats

  rx_bits = []; rx_np = []; timing_est_log = []; delta_t_log = []; foff_est_hz_log = [];
  phase_est_pilot_log = [];
  Nerrs = Nbits = 0;

  % load real samples from file

  Ascale= 4E5;
  frx=fopen(filename,"rb"); rx = 2*fread(frx, Inf, "short")/4E5; fclose(frx);
  Nsam = length(rx); Nframes = floor(Nsam/Nsamperframe);
  prx = 1;

  % 'prime' rx buf to get correct coarse timing (for now)

  prx = 1;
  states.rxbuf(M+Ncp+2*Nsamperframe+1:Nrxbuf) = rx(prx:Nsamperframe+2*(M+Ncp));
  prx += Nsamperframe+2*(M+Ncp);

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

    [rx_bits states aphase_est_pilot_log arx_np] = ofdm_demod(states, rxbuf_in);

    rx_np = [rx_np arx_np];
    timing_est_log = [timing_est_log states.timing_est];
    delta_t_log = [delta_t_log states.delta_t];
    foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
    phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];

    % measure bit errors

    errors = xor(tx_bits, rx_bits);
    Nerrs += sum(errors);
    Nerrs_log(f) = sum(errors);
    Nbits += Nbitsperframe;
  end

  printf("BER: %5.4f Nbits: %d Nerrs: %d\n", Nerrs/Nbits, Nbits, Nerrs);

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
  axis([1 max(Nframes,2) -3 3]);
  title('Fine Freq');

  figure(5); clf;
  %plot(Nerrs_log);
  
endfunction
