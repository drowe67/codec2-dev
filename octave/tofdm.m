% tofdm.m
% David Rowe and Steve Sampson June 2017
%
% Octave script for comparing Octave and C versions of OFDZM modem

% ------------------------------------------------------------------
1;

function rx = ebno_awgn_channel(tx, Fs, Rb, EbNodB)
  EbNo = 10^(EbNodB/10);
  variance = Fs / (Rb*EbNo);
  noise = sqrt(variance)*randn(1,length(tx));
  avg = sum(abs(tx))/length(tx);
  rx = noise.*avg + tx;
  %rx = rx * .10;
end

function nums = im_re_interleave(nim)
  nums = zeros(1,length(nim)*2);
  nums(1:2:length(nums)) = real(nim);
  nums(2:2:length(nums)) = imag(nim);

end

function pass = run_ofdm_test(Nframes,sample_clock_offset_ppm,foff_hz,EbNodB = 100,nheadsamp = 0)
  %Nframes = 30;
  %sample_clock_offset_ppm = 100;
  %foff_hz = 5;

  more off; format;
  ofdm_lib;
  autotest;

  % ---------------------------------------------------------------------
  % Run Octave version 
  % ---------------------------------------------------------------------

  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 16; Ns = 8;
  Rb_real = 700;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;

  rand('seed',2);
  tx_bits = round(rand(1,Nbitsperframe));

  % Run tx loop

  tx_bits_log = []; tx_log = [];
  for f=1:Nframes
    tx_bits_log = [tx_bits_log tx_bits];
    tx_log = [tx_log ofdm_mod(states, tx_bits)];
  end

  % Channel simulation ----------------------------------------------

  rx_log = sample_clock_offset(tx_log, sample_clock_offset_ppm);
  rx_log = freq_shift(rx_log, foff_hz, Fs);
  #rx_log(1:nheadsamp) = zeros(nheadsamp,1);
  rx_log_l = length(rx_log)
  rx_log = [zeros(1,nheadsamp), rx_log];
  rx_log = rx_log(1:rx_log_l);
  rx_log = ebno_awgn_channel(rx_log, Fs, Rb_real, EbNodB);

  rx_vec = fopen("tofdm_rx_vec","wb+");
  fwrite(rx_vec,im_re_interleave(rx_log),"single");
  fclose(rx_vec);

  % Rx ---------------------------------------------------------------

  % Init rx with ideal timing so we can test with timing estimation disabled

  Nsam = length(rx_log);
  prx = 1;
  nin = Nsamperframe+2*(M+Ncp);
  states.rxbuf(Nrxbuf-nin+1:Nrxbuf) = rx_log(prx:nin);
  prx += nin;

  rxbuf_log = []; rxbuf_in_log = []; rx_sym_log = []; foff_hz_log = []; 
  timing_est_log = []; sample_point_log = []; 
  phase_est_pilot_log = []; rx_amp_log = [];
  rx_np_log = []; rx_bits_log = [];

  states.timing_en = 1;
  states.foff_est_en = 1;
  states.phase_est_en = 1;

  if states.timing_en == 0
    % manually set ideal timing instant
    states.sample_point = Ncp;
  end

  %states.foff_est_hz = 10;
  for f=1:Nframes

    % insert samples at end of buffer, set to zero if no samples
    % available to disable phase estimation on future pilots on last
    % frame of simulation
  
    nin = states.nin;
    lnew = min(Nsam-prx+1,nin);
    rxbuf_in = zeros(1,nin);
    %printf("nin: %d prx: %d lnew: %d\n", nin, prx, lnew);
    if lnew
      rxbuf_in(1:lnew) = rx_log(prx:prx+lnew-1);
    end
    prx += lnew;

    [rx_bits states aphase_est_pilot_log arx_np arx_amp] = ofdm_demod2(states, rxbuf_in);

    % log some states for comparison to C

    rxbuf_in_log = [rxbuf_in_log rxbuf_in];
    rxbuf_log = [rxbuf_log states.rxbuf];
    rx_sym_log = [rx_sym_log; states.rx_sym];
    phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];
    rx_amp_log = [rx_amp_log arx_amp];
    foff_hz_log = [foff_hz_log; states.foff_est_hz];
    timing_est_log = [timing_est_log; states.timing_est];
    sample_point_log = [sample_point_log; states.sample_point];
    rx_np_log = [rx_np_log arx_np];
    rx_bits_log = [rx_bits_log rx_bits];
    
  end

  % ---------------------------------------------------------------------
  % Run C version and plot Octave and C states and differences 
  % ---------------------------------------------------------------------

  % Override default path by setting path_to_tofdm = "/your/path/to/tofdm"

  if exist("path_to_tofdm", "var") == 0
    path_to_tofdm = "../build_linux/unittest/tofdm";
  end
  system(path_to_tofdm);

  load tofdm_out.txt;
  % Generated with modem probe thing
  load ofdm_test.txt;

  fg = 1;
  figure(fg++); clf; plot(rx_np_log,'+'); title('Octave Scatter Diagram'); axis([-1.5 1.5 -1.5 1.5]);
  figure(fg++); clf; plot(rx_np_log_c,'+'); title('C Scatter Diagram'); axis([-1.5 1.5 -1.5 1.5]);

  stem_sig_and_error(fg++, 111, tx_bits_log_c, tx_bits_log - tx_bits_log_c, 'tx bits', [1 length(tx_bits_log) -1.5 1.5])

  stem_sig_and_error(fg, 211, real(tx_log_c), real(tx_log - tx_log_c), 'tx re', [1 length(tx_log_c) -0.1 0.1])
  stem_sig_and_error(fg++, 212, imag(tx_log_c), imag(tx_log - tx_log_c), 'tx im', [1 length(tx_log_c) -0.1 0.1])

  stem_sig_and_error(fg, 211, real(rx_log_c), real(rx_log - rx_log_c), 'rx re', [1 length(rx_log_c) -0.1 0.1])
  stem_sig_and_error(fg++, 212, imag(rx_log_c), imag(rx_log - rx_log_c), 'rx im', [1 length(rx_log_c) -0.1 0.1])

  stem_sig_and_error(fg, 211, real(rxbuf_in_log_c), real(rxbuf_in_log - rxbuf_in_log_c), 'rxbuf in re', [1 length(rxbuf_in_log_c) -0.1 0.1])
  stem_sig_and_error(fg++, 212, imag(rxbuf_in_log_c), imag(rxbuf_in_log - rxbuf_in_log_c), 'rxbuf in im', [1 length(rxbuf_in_log_c) -0.1 0.1])

  stem_sig_and_error(fg, 211, real(rxbuf_log_c), real(rxbuf_log - rxbuf_log_c), 'rxbuf re', [1 length(rxbuf_log_c) -0.1 0.1])
  stem_sig_and_error(fg++, 212, imag(rxbuf_log_c), imag(rxbuf_log - rxbuf_log_c), 'rxbuf im', [1 length(rxbuf_log_c) -0.1 0.1])

  stem_sig_and_error(fg, 211, real(rx_sym_log_c), real(rx_sym_log - rx_sym_log_c), 'rx sym re', [1 length(rx_sym_log_c) -1.5 1.5])
  stem_sig_and_error(fg++, 212, imag(rx_sym_log_c), imag(rx_sym_log - rx_sym_log_c), 'rx sym im', [1 length(rx_sym_log_c) -1.5 1.5])

  % for angles pi and -pi are the same

  d = phase_est_pilot_log - phase_est_pilot_log_c; d = angle(exp(j*d));

  stem_sig_and_error(fg, 211, phase_est_pilot_log_c, d, 'phase est pilot', [1 length(phase_est_pilot_log_c) -1.5 1.5])
  stem_sig_and_error(fg++, 212, rx_amp_log_c, rx_amp_log - rx_amp_log_c, 'rx amp', [1 length(rx_amp_log_c) -1.5 1.5])

  stem_sig_and_error(fg++, 111, foff_hz_log_c, (foff_hz_log - foff_hz_log_c), 'foff hz', [1 length(foff_hz_log_c) -1.5 1.5])

  stem_sig_and_error(fg,   211, timing_est_log_c, (timing_est_log - timing_est_log_c), 'timing est', [1 length(timing_est_log_c) -1.5 1.5])
  stem_sig_and_error(fg++, 212, sample_point_log_c, (sample_point_log - sample_point_log_c), 'sample point', [1 length(sample_point_log_c) -1.5 1.5])

  stem_sig_and_error(fg++, 111, rx_bits_log_c, rx_bits_log - rx_bits_log_c, 'rx bits', [1 length(rx_bits_log) -1.5 1.5])

  % Run through checklist -----------------------------
  pass = true;
  pass = check_no_abs(W, W_c, 'W') && pass;
  pass = check(tx_bits_log, tx_bits_log_c, 'tx_bits') && pass;
  pass = check(tx_log, tx_log_c, 'tx') && pass;
  pass = check(rx_log, rx_log_c, 'rx') && pass;
  pass = check(rxbuf_in_log, rxbuf_in_log_c, 'rxbuf in') && pass;
  pass = check(rxbuf_log, rxbuf_log_c, 'rxbuf') && pass;
  pass = check(rx_sym_log, rx_sym_log_c, 'rx_sym') && pass;
  pass = check(phase_est_pilot_log, phase_est_pilot_log_c, 'phase_est_pilot', tol=1E-3, its_an_angle=1) && pass;
  pass = check(rx_amp_log, rx_amp_log_c, 'rx_amp') && pass;
  pass = check(timing_est_log, timing_est_log_c', 'timing_est') && pass;
  pass = check(sample_point_log, sample_point_log_c, 'sample_point') && pass;
  pass = check(foff_hz_log, foff_hz_log_c', 'foff_est_hz') && pass;
  pass = check(rx_bits_log, rx_bits_log_c, 'rx_bits') && pass;
  
  figure(fg++)
  stem(timing_est_log-timing_est_log_c')
  figure(fg++)
  stem(phase_est_pilot_log-phase_est_pilot_log_c)

end

% This works best on my machine -- Brady
graphics_toolkit('fltk')

run_ofdm_test(60,20,0.5,40,100)