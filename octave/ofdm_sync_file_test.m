% ofdm_sync_file_test.m
% David Rowe April 2017
%
% Reads an off air file, and dumps sync metrics.  Used for sync debug,
% also see other functions on ofdm_dev

function ofdm_sync_file_test(filename, mode="700D")
  ofdm_lib;
  more off;

  % init modem

  [bps Rs Tcp Ns Nc] = ofdm_init_mode(mode);
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;
  states.verbose = 1;

  % load real samples from file

  Ascale= states.amp_scale/2.0;  % /2 as real signal has half amplitude
  frx=fopen(filename,"rb"); rx = fread(frx, Inf, "short")/Ascale; fclose(frx);
  Nsam = length(rx); Nframes = floor(Nsam/Nsamperframe);
  prx = 1;

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
    
    states.rxbuf(1:Nrxbuf-states.nin) = states.rxbuf(states.nin+1:Nrxbuf);
    states.rxbuf(Nrxbuf-states.nin+1:Nrxbuf) = rxbuf_in;
    st = M+Ncp + Nsamperframe + 1; en = st + 2*Nsamperframe;
    #{
    size(states.rxbuf(st:en))
    [ct_est timing_valid timing_mx] = est_timing(states, states.rxbuf(st:en), states.rate_fs_pilot_samples);
    [foff_est states] = est_freq_offset(states, states.rxbuf(st:en), states.rate_fs_pilot_samples, ct_est);
    #}
    st = (f-1)*Nsamperframe+1; en = st + 2*Nsamperframe;
    [ct_est timing_valid timing_mx] = est_timing(states, rx(st:en)', states.rate_fs_pilot_samples);
    [foff_est states] = est_freq_offset(states, rx(st:en)', states.rate_fs_pilot_samples, ct_est);
    
    printf("  ct_est: %d mx: %3.2f coarse_foff: %4.1f\n", ct_est, timing_mx, foff_est);

    %[timing_valid states] = ofdm_sync_search(states, rxbuf_in);
  end

endfunction
