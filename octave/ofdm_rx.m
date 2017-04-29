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

  % load real samples from file

  frx=fopen(filename,"rb"); rx = fread(frx, Inf, "short"); fclose(frx);

  Nsam = length(rx); Nframes = floor(Nsam/Nsamperframe);
  prx = 1;
  rx_bits = []; rx_np = []; timing_est_log = []; delta_t_log = []; foff_est_hz_log = [];
  phase_est_pilot_log = [];

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

    [arx_bits states aphase_est_pilot_log arx_np] = ofdm_demod(states, rxbuf_in);

    rx_bits = [rx_bits arx_bits]; rx_np = [rx_np arx_np];
    timing_est_log = [timing_est_log states.timing_est];
    delta_t_log = [delta_t_log states.delta_t];
    foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
    phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];
  end

  figure(1); clf; 
  plot(rx_np,'+');
  %axis([-2 2 -2 2]);
  title('Scatter');

endfunction
