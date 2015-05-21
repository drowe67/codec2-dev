% tcohpsk.m
% David Rowe Oct 2014
%
% Octave script that:
%
% i) tests the C port of the coherent PSK modem.  This script loads
%    the output of unittest/tcohpsk.c and compares it to the output of
%    the reference versions of the same functions written in Octave.
%
% or (ii) can optionally be used to run an Octave version of the cohpsk
%    modem to tune and develop it.

%  TODO:
%
%  [ ] Test
%      [X] AWGN channel
%      [X] freq offset
%      [X] fading channel
%      [X] freq drift
%      [ ] timing drift
%  [X] tune perf/impl loss to get closer to ideal
%      [X] linear interp of phase for better fading perf
%  [X] freq offset/drift feedback loop 
%  [X] PAPR measurement and reduction
%  [X] false sync
%      [X] doesn't sync up on noise (used EsNo = -12)
%      [X] similar but invalid signal like huge f off
%  [X] ability to "unsync" when signal dissapears
%  [ ] some calibrated tests against FreeDV 1600
%      + compare sound quality at various Es/Nos
%  [ ] way to handle eom w/o nasties
%      + like mute ouput when signal has gone or v low snr
%      + instantaneous snr
%  [ ] nasty rig filter passband

graphics_toolkit ("gnuplot");
more off;

cohpsk;
fdmdv;
autotest;

rand('state',1); 
randn('state',1);

% select which test  ----------------------------------------------------------

%test = 'compare to c';
test = 'awgn';
%test = 'fading';

% some parameters that can be over ridden, e.g. to disable parts of modem

initial_sync = 0;  % setting this to 1 put us straight into sync w/o freq offset est
ftrack_en    = 1;  % set to 0 to disable freq tracking

% predefined tests ....

if strcmp(test, 'compare to c')
  frames = 35;
  foff =  55.5;
  dfoff = -0.5/Fs;
  EsNodB = 8;
  fading_en = 0;
  hf_delay_ms = 2;
  compare_with_c = 1;
end

% should be BER around 0.015 to 0.02

if strcmp(test, 'awgn')
  frames = 10;
  foff =  0;
  dfoff = 0;
  EsNodB = 80;
  fading_en = 0;
  hf_delay_ms = 2;
  compare_with_c = 0;
end

% Similar to AWGN - should be BER around 0.015 to 0.02

if strcmp(test, 'fading');
  frames = 100;
  foff = -52.3;
  dfoff = 0.5/Fs;
  EsNodB = 12;
  fading_en = 1;
  hf_delay_ms = 2;
  compare_with_c = 0;
end

EsNo = 10^(EsNodB/10);

% modem constants ----------------------------------------------------------

Rs = 75;               % symbol rate in Hz
Nc = 4;                % number of carriers
Nd = 2;                % diveristy factor
framesize = 32;        % number of payoad data bits in the frame

Nsw = 3;               % frames we demod for initial sync window
afdmdv.Nsym = 6;       % size of tx/tx root nyquist filter in symbols
afdmdv.Nt = 5;         % number of symbols we estimate timing over

% FDMDV init ---------------------------------------------------------------

Fs = afdmdv.Fs = 7500;
afdmdv.Nc = Nd*Nc-1;
afdmdv.Rs = Rs;
M = afdmdv.M  = afdmdv.Fs/afdmdv.Rs;
afdmdv.Nfilter = afdmdv.Nsym*M;
afdmdv.tx_filter_memory = zeros(afdmdv.Nc+1, afdmdv.Nfilter);
excess_bw = 0.5;
afdmdv.gt_alpha5_root = gen_rn_coeffs(excess_bw, 1/Fs, Rs, afdmdv.Nsym, afdmdv.M);

Fcentre = afdmdv.Fcentre = 1500;
afdmdv.Fsep = afdmdv.Rs*(1+excess_bw);
afdmdv.phase_tx = ones(afdmdv.Nc+1,1);
freq_hz = afdmdv.Fsep*( -Nc*Nd/2 - 0.5 + (1:Nc*Nd).^0.98 );
afdmdv.freq_pol = 2*pi*freq_hz/Fs;
afdmdv.freq = exp(j*afdmdv.freq_pol);
afdmdv.Fcentre = 1500;

afdmdv.fbb_rect = exp(j*2*pi*Fcentre/Fs);
afdmdv.fbb_phase_tx = 1;
afdmdv.fbb_phase_rx = 1;

afdmdv.Nrxdec = 31;
afdmdv.rxdec_coeff = fir1(afdmdv.Nrxdec-1, 0.25)';
afdmdv.rxdec_lpf_mem = zeros(1,afdmdv.Nrxdec-1+afdmdv.M);

P = afdmdv.P = 4;
afdmdv.phase_rx = ones(afdmdv.Nc+1,1);
afdmdv.Nfilter = afdmdv.Nsym*afdmdv.M;
afdmdv.rx_fdm_mem = zeros(1,afdmdv.Nfilter + afdmdv.M);
Q = afdmdv.Q = afdmdv.M/4;

afdmdv.rx_filter_mem_timing = zeros(afdmdv.Nc+1, afdmdv.Nt*afdmdv.P);
afdmdv.Nfiltertiming = afdmdv.M + afdmdv.Nfilter + afdmdv.M;

afdmdv.rx_filter_memory = zeros(afdmdv.Nc+1, afdmdv.Nfilter);

afdmdv.filt = 0;
afdmdv.prev_rx_symb = ones(1,afdmdv.Nc+1);

% COHPSK Init --------------------------------------------------------

acohpsk = standard_init();
acohpsk.framesize        = framesize;
acohpsk.ldpc_code        = 0;
acohpsk.ldpc_code_rate   = 1;
acohpsk.Nc               = Nc;
acohpsk.Rs               = Rs;
acohpsk.Ns               = 4;
acohpsk.coh_en           = 1;
acohpsk.Nd               = Nd;
acohpsk.modulation       = 'qpsk';
acohpsk.do_write_pilot_file = 0;
acohpsk = symbol_rate_init(acohpsk);
acohpsk.Ndft = 1024;
acohpsk.f_est = afdmdv.Fcentre;

ch_fdm_frame_buf = zeros(1, Nsw*acohpsk.Nsymbrowpilot*afdmdv.M);

% -----------------------------------------------------------

tx_bits_log = [];
tx_symb_log = [];
rx_amp_log = [];
rx_phi_log = [];
ch_symb_log = [];
rx_symb_log = [];
rx_bits_log = [];
tx_bits_prev_log = [];
uvnoise_log = [];
nerr_log = [];
tx_baseband_log = [];
tx_fdm_frame_log = [];
ch_fdm_frame_log = [];
rx_fdm_frame_bb_log = [];
rx_filt_log = [];
rx_fdm_filter_log = [];
rx_baseband_log = [];
rx_fdm_frame_log = [];
ct_symb_ff_log = [];
rx_timing_log = [];
ratio_log = [];
foff_log = [];
f_est_log = [];

% Channel modeling and BER measurement ----------------------------------------

rand('state',1); 
tx_bits_coh = round(rand(1,framesize*10));
ptx_bits_coh = 1;

Nerrs = Tbits = 0;
prev_tx_bits = [];

phase_ch = 1;
sync = initial_sync;
acohpsk.f_est = Fcentre;
acohpsk.f_fine_est = 0;
acohpsk.ct = 4;
acohpsk.ftrack_en = ftrack_en;

[spread spread_2ms hf_gain] = init_hf_model(Fs, Fs, frames*acohpsk.Nsymbrowpilot*afdmdv.M);
hf_n = 1;
nhfdelay = floor(hf_delay_ms*Fs/1000);
ch_fdm_delay = zeros(1, acohpsk.Nsymbrowpilot*M + nhfdelay);

% main loop --------------------------------------------------------------------

for f=1:frames
  acohpsk.frame = f;
  tx_bits = tx_bits_coh(ptx_bits_coh:ptx_bits_coh+framesize-1);
  ptx_bits_coh += framesize;
  if ptx_bits_coh > length(tx_bits_coh)
    ptx_bits_coh = 1;
  end

  tx_bits_log = [tx_bits_log tx_bits];

  [tx_symb tx_bits] = bits_to_qpsk_symbols(acohpsk, tx_bits, [], []);
  tx_symb_log = [tx_symb_log; tx_symb];

  tx_fdm_frame = [];
  for r=1:acohpsk.Nsymbrowpilot
    tx_onesymb = tx_symb(r,:);
    [tx_baseband afdmdv] = tx_filter(afdmdv, tx_onesymb);
    tx_baseband_log = [tx_baseband_log tx_baseband];
    [tx_fdm afdmdv] = fdm_upconvert(afdmdv, tx_baseband);
    tx_fdm_frame = [tx_fdm_frame tx_fdm];
  end

  clip = 5;
  ind = find(abs(tx_fdm_frame) > clip);
  tx_fdm_frame(ind) = clip*exp(j*angle(tx_fdm_frame(ind)));

  tx_fdm_frame_log = [tx_fdm_frame_log tx_fdm_frame];

  %
  % Channel --------------------------------------------------------------------
  %

  %[ch_fdm_frame phase_ch] = freq_shift(tx_fdm_frame, foff, Fs, phase_ch);
  %foff += dfoff*acohpsk.Nsymbrowpilot*M;

  ch_fdm_frame = zeros(1,acohpsk.Nsymbrowpilot*M);
  for i=1:acohpsk.Nsymbrowpilot*M
    foff_rect = exp(j*2*pi*foff/Fs);
    foff += dfoff;
    phase_ch *= foff_rect;
    ch_fdm_frame(i) = tx_fdm_frame(i) * phase_ch;
  end
  foff_log = [foff_log foff];
  phase_ch /= abs(phase_ch);
  % printf("foff: %f  ", foff);

  if fading_en
    ch_fdm_delay(1:nhfdelay) = ch_fdm_delay(acohpsk.Nsymbrowpilot*M+1:nhfdelay+acohpsk.Nsymbrowpilot*M);
    ch_fdm_delay(nhfdelay+1:nhfdelay+acohpsk.Nsymbrowpilot*M) = ch_fdm_frame;

    for i=1:acohpsk.Nsymbrowpilot*M
      ahf_model = hf_gain*(spread(hf_n)*ch_fdm_frame(i) + spread_2ms(hf_n)*ch_fdm_delay(i));
      ch_fdm_frame(i) = ahf_model;
      hf_n++;
    end
  end

  % each carrier has power = 2, total power 2Nc, total symbol rate NcRs, noise BW B=Fs
  % Es/No = (C/Rs)/(N/B), N = var = 2NcFs/NcRs(Es/No) = 2Fs/Rs(Es/No)

  variance = 2*Fs/(acohpsk.Rs*EsNo);
  uvnoise = sqrt(0.5)*(randn(1,acohpsk.Nsymbrowpilot*M) + j*randn(1,acohpsk.Nsymbrowpilot*M));
  uvnoise_log = [uvnoise_log uvnoise];
  noise = sqrt(variance)*uvnoise;

  ch_fdm_frame += noise;

  ch_fdm_frame_log = [ch_fdm_frame_log ch_fdm_frame];

  %
  % Demod ----------------------------------------------------------------------
  %

  % store two frames of received samples so we can rewind if we get a good candidate

  ch_fdm_frame_buf(1:(Nsw-1)*acohpsk.Nsymbrowpilot*M) = ch_fdm_frame_buf(acohpsk.Nsymbrowpilot*M+1:Nsw*acohpsk.Nsymbrowpilot*M);
  ch_fdm_frame_buf((Nsw-1)*acohpsk.Nsymbrowpilot*M+1:Nsw*acohpsk.Nsymbrowpilot*M) = ch_fdm_frame;

  next_sync = sync;
  nin = M;

  % if out of sync do Initial Freq offset estimation over NSW frames to flush out memories

  if (sync == 0)

    % we can test +/- 20Hz, so we break this up into 3 tests to cover +/- 60Hz

    max_ratio = 0;
    for acohpsk.f_est = Fcentre-40:40:Fcentre+40
        
      printf("  [%d] acohpsk.f_est: %f +/- 20\n", f, acohpsk.f_est);

      % we are out of sync so reset f_est and process two frames to clean out memories

      [ch_symb rx_timing rx_filt rx_baseband afdmdv acohpsk.f_est] = rate_Fs_rx_processing(afdmdv, ch_fdm_frame_buf, acohpsk.f_est, Nsw*acohpsk.Nsymbrowpilot, nin, 0);
      rx_baseband_log = [rx_baseband_log rx_baseband];
 
      rx_filt_log = [rx_filt_log rx_filt];
      ch_symb_log = [ch_symb_log; ch_symb];

      for i=1:Nsw-1
        acohpsk.ct_symb_buf = update_ct_symb_buf(acohpsk.ct_symb_buf, ch_symb((i-1)*acohpsk.Nsymbrowpilot+1:i*acohpsk.Nsymbrowpilot,:), acohpsk.Nct_sym_buf, acohpsk.Nsymbrowpilot);
      end
      [anext_sync acohpsk] = frame_sync_fine_freq_est(acohpsk, ch_symb((Nsw-1)*acohpsk.Nsymbrowpilot+1:Nsw*acohpsk.Nsymbrowpilot,:), sync, next_sync);

      if anext_sync == 1
        %printf("  [%d] acohpsk.ratio: %f\n", f, acohpsk.ratio);
        if acohpsk.ratio > max_ratio
          max_ratio   = acohpsk.ratio;
          f_est       = acohpsk.f_est - acohpsk.f_fine_est;
          next_sync   = anext_sync;
        end
      end
    end

    if next_sync == 1

      % we've found a sync candidate!
      % re-process last two frames with adjusted f_est then check again

      acohpsk.f_est = f_est;

      printf("  [%d] trying sync and f_est: %f\n", f, acohpsk.f_est);

      [ch_symb rx_timing rx_filt rx_baseband afdmdv f_est] = rate_Fs_rx_processing(afdmdv, ch_fdm_frame_buf, acohpsk.f_est, Nsw*acohpsk.Nsymbrowpilot, nin, 0);
      rx_baseband_log = [rx_baseband_log rx_baseband];
      rx_filt_log = [rx_filt_log rx_filt];
      ch_symb_log = [ch_symb_log; ch_symb];
 
      for i=1:Nsw-1
        acohpsk.ct_symb_buf = update_ct_symb_buf(acohpsk.ct_symb_buf, ch_symb((i-1)*acohpsk.Nsymbrowpilot+1:i*acohpsk.Nsymbrowpilot,:), acohpsk.Nct_sym_buf, acohpsk.Nsymbrowpilot);
      end
      [next_sync acohpsk] = frame_sync_fine_freq_est(acohpsk, ch_symb((Nsw-1)*acohpsk.Nsymbrowpilot+1:Nsw*acohpsk.Nsymbrowpilot,:), sync, next_sync);
      if abs(acohpsk.f_fine_est) > 2
        printf("  [%d] Hmm %f is a bit big so back to coarse est ...\n", f, acohpsk.f_fine_est);
        next_sync = 0;
      end

      if acohpsk.ratio < 0.9
        next_sync = 0;
      end
      if next_sync == 1
        % OK we are in sync!
        % demodulate first frame (demod completed below)

        printf("  [%d] in sync! f_est: %f ratio: %f \n", f, f_est, acohpsk.ratio);
        acohpsk.ct_symb_ff_buf(1:acohpsk.Nsymbrowpilot+2,:) = acohpsk.ct_symb_buf(acohpsk.ct+1:acohpsk.ct+acohpsk.Nsymbrowpilot+2,:);
      end
    end  
  end

  % If in sync just do sample rate processing on latest frame

  if sync == 1
    [ch_symb rx_timing rx_filt rx_baseband afdmdv acohpsk.f_est] = rate_Fs_rx_processing(afdmdv, ch_fdm_frame, acohpsk.f_est, acohpsk.Nsymbrowpilot, nin, acohpsk.ftrack_en);
    [next_sync acohpsk] = frame_sync_fine_freq_est(acohpsk, ch_symb, sync, next_sync);

    acohpsk.ct_symb_ff_buf(1:2,:) = acohpsk.ct_symb_ff_buf(acohpsk.Nsymbrowpilot+1:acohpsk.Nsymbrowpilot+2,:);
    acohpsk.ct_symb_ff_buf(3:acohpsk.Nsymbrowpilot+2,:) = acohpsk.ct_symb_buf(acohpsk.ct+3:acohpsk.ct+acohpsk.Nsymbrowpilot+2,:);

    rx_baseband_log = [rx_baseband_log rx_baseband];
    rx_filt_log = [rx_filt_log rx_filt];
    ch_symb_log = [ch_symb_log; ch_symb];     
    rx_timing_log = [rx_timing_log rx_timing];
    f_est_log = [f_est_log acohpsk.f_est];
    %printf("%f\n", acohpsk.f_est);
  end

  % if we are in sync complete demodulation with symbol rate processing

  if (next_sync == 1) || (sync == 1)
    [rx_symb rx_bits rx_symb_linear amp_ phi_ EsNo_] = qpsk_symbols_to_bits(acohpsk, acohpsk.ct_symb_ff_buf);
    rx_symb_log = [rx_symb_log; rx_symb];
    rx_amp_log = [rx_amp_log; amp_];
    rx_phi_log = [rx_phi_log; phi_];
    rx_bits_log = [rx_bits_log rx_bits];
    tx_bits_prev_log = [tx_bits_prev_log prev_tx_bits2];
    ratio_log = [ratio_log acohpsk.ratio];
    ct_symb_ff_log = [ct_symb_ff_log; acohpsk.ct_symb_ff_buf(1:acohpsk.Nsymbrowpilot,:)];

    % BER stats

    if f > 2
      error_positions = xor(prev_tx_bits2, rx_bits);
      %error_positions = xor(prev_tx_bits, rx_bits);
      Nerrs  += sum(error_positions);
      nerr_log = [nerr_log sum(error_positions)];
      Tbits += length(error_positions);
    end
    printf("\r  [%d]", f);
  end

   
  %rx_fdm_frame_bb_log = [rx_fdm_frame_bb_log rx_fdm_frame_bb];
  %rx_baseband_log = [rx_baseband_log rx_baseband];
  %rx_filt_log = [rx_filt_log rx_filt];
  %rx_timing_log = [rx_timing_log rx_timing];
  %ch_symb_log = [ch_symb_log; ch_symb];
  % ct_symb_ff_log = [ct_symb_ff_log; acohpsk.ct_symb_ff_buf(1:acohpsk.Nsymbrowpilot,:)];


  if sync == 0
    Nerrs = 0;
    Tbits = 0;
    nerr_log = [];
  end

  % printf("f: %d sync: %d next_sync: %d\n", f, sync, next_sync);
  [sync acohpsk] = sync_state_machine(acohpsk, sync, next_sync);

  prev_tx_bits2 = prev_tx_bits;
  prev_tx_bits = tx_bits;

end

ber = Nerrs/Tbits;
printf("\nOctave EsNodB: %4.1f ber..: %4.3f Nerrs..: %d Tbits..: %d\n", EsNodB, ber, Nerrs, Tbits);

if compare_with_c

  % Output vectors from C port ---------------------------------------------------

  load ../build_linux/unittest/tcohpsk_out.txt

  % Determine bit error rate

  sz = length(tx_bits_log_c);
  Nerrs_c = sum(xor(tx_bits_prev_log, rx_bits_log_c(framesize+1:length(rx_bits_log_c))));
  Tbits_c = length(tx_bits_prev_log);
  ber_c = Nerrs_c/Tbits_c;
  printf("C EsNodB.....: %4.1f ber_c: %4.3f Nerrs_c: %d Tbits_c: %d\n", EsNodB, ber_c, Nerrs_c, Tbits_c);

  stem_sig_and_error(1, 111, tx_bits_log_c, tx_bits_log - tx_bits_log_c, 'tx bits', [1 length(tx_bits) -1.5 1.5])
  stem_sig_and_error(2, 211, real(tx_symb_log_c), real(tx_symb_log - tx_symb_log_c), 'tx symb re', [1 length(tx_symb_log_c) -1.5 1.5])
  stem_sig_and_error(2, 212, imag(tx_symb_log_c), imag(tx_symb_log - tx_symb_log_c), 'tx symb im', [1 length(tx_symb_log_c) -1.5 1.5])

  stem_sig_and_error(3, 211, real(tx_fdm_frame_log_c), real(tx_fdm_frame_log - tx_fdm_frame_log_c), 'tx fdm frame re', [1 length(tx_fdm_frame_log) -10 10])
  stem_sig_and_error(3, 212, imag(tx_fdm_frame_log_c), imag(tx_fdm_frame_log - tx_fdm_frame_log_c), 'tx fdm frame im', [1 length(tx_fdm_frame_log) -10 10])
  stem_sig_and_error(4, 211, real(ch_fdm_frame_log_c), real(ch_fdm_frame_log - ch_fdm_frame_log_c), 'ch fdm frame re', [1 length(ch_fdm_frame_log) -10 10])
  stem_sig_and_error(4, 212, imag(ch_fdm_frame_log_c), imag(ch_fdm_frame_log - ch_fdm_frame_log_c), 'ch fdm frame im', [1 length(ch_fdm_frame_log) -10 10])

  stem_sig_and_error(5, 211, real(rx_baseband_log_c(1,:)), real(rx_baseband_log(1,:) - rx_baseband_log_c(1,:)), 'rx baseband re', [1 length(rx_baseband_log) -10 10])
  stem_sig_and_error(5, 212, imag(rx_baseband_log_c(1,:)), imag(rx_baseband_log(1,:) - rx_baseband_log_c(1,:)), 'rx baseband im', [1 length(rx_baseband_log) -10 10])

  [n m] = size(ch_symb_log);
  stem_sig_and_error(6, 211, real(ch_symb_log_c), real(ch_symb_log - ch_symb_log_c), 'ch symb re', [1 n -1.5 1.5])
  stem_sig_and_error(6, 212, imag(ch_symb_log_c), imag(ch_symb_log - ch_symb_log_c), 'ch symb im', [1 n -1.5 1.5])

  [n m] = size(rx_symb_log);
  stem_sig_and_error(8, 211, rx_amp_log_c, rx_amp_log - rx_amp_log_c, 'Amp Est', [1 n -1.5 1.5])
  stem_sig_and_error(8, 212, rx_phi_log_c, rx_phi_log - rx_phi_log_c, 'Phase Est', [1 n -4 4])
  stem_sig_and_error(9, 211, real(rx_symb_log_c), real(rx_symb_log - rx_symb_log_c), 'rx symb re', [1 n -1.5 1.5])
  stem_sig_and_error(9, 212, imag(rx_symb_log_c), imag(rx_symb_log - rx_symb_log_c), 'rx symb im', [1 n -1.5 1.5])

  stem_sig_and_error(10, 111, rx_bits_log_c, rx_bits_log - rx_bits_log_c, 'rx bits', [1 length(rx_bits_log) -1.5 1.5])
  stem_sig_and_error(11, 111, f_est_log_c - Fcentre, f_est_log - f_est_log_c, 'f est', [1 length(f_est_log) foff-5 foff+5])

  check(tx_bits_log, tx_bits_log_c, 'tx_bits');
  check(tx_symb_log, tx_symb_log_c, 'tx_symb');
  check(tx_fdm_frame_log, tx_fdm_frame_log_c, 'tx_fdm_frame');
  check(ch_fdm_frame_log, ch_fdm_frame_log_c, 'ch_fdm_frame');
  %check(rx_fdm_frame_bb_log, rx_fdm_frame_bb_log_c, 'rx_fdm_frame_bb', 0.01);

  check(ch_symb_log, ch_symb_log_c, 'ch_symb',0.05);
  %check(ct_symb_ff_log, ct_symb_ff_log_c, 'ct_symb_ff',0.01);
  check(rx_amp_log, rx_amp_log_c, 'rx_amp_log',0.01);
  check(rx_phi_log, rx_phi_log_c, 'rx_phi_log',0.05);
  check(rx_symb_log, rx_symb_log_c, 'rx_symb',0.01);
  check(rx_bits_log, rx_bits_log_c, 'rx_bits');
  check(f_est_log, f_est_log_c, 'f_est');


else
  
  papr = max(tx_fdm_frame_log.*conj(tx_fdm_frame_log)) / mean(tx_fdm_frame_log.*conj(tx_fdm_frame_log));
  papr_dB = 10*log10(papr);
  printf("av tx pwr: %f PAPR: %4.2f dB av rx pwr: %f\n", var(tx_fdm_frame_log), papr_dB, var(ch_fdm_frame_log));

  % some other useful plots

  figure(1)
  clf
  subplot(211)
  plot(real(tx_fdm_frame_log))
  title('tx fdm real');
  subplot(212)
  plot(imag(tx_fdm_frame_log))
  title('tx fdm imag');

  figure(2)
  clf;
  plot(rx_symb_log*exp(j*pi/4),'+')
  title('Scatter');

  figure(3)
  clf;
  plot(rx_timing_log)
  title('rx timing');

  figure(4)
  clf;
  stem(nerr_log)
  title('Bit Errors');

  figure(5)
  clf;
  stem(ratio_log)
  title('Sync ratio');

  figure(6);
  clf;
  subplot(211)
  plot(foff_log,';freq offset;');
  hold on;
  plot(f_est_log - Fcentre,'g;freq offset est;');
  hold off;
  title('freq offset');
  legend("boxoff");  
  subplot(212)
  plot(foff_log(1:length(f_est_log)) - f_est_log + Fcentre)
  title('freq offset estimation error');

  figure(7)
  clf
  spec = 20*log10(abs(fft(tx_fdm_frame_log)));
  l = length(spec);
  plot((Fs/l)*(1:l), spec)
  axis([1 Fs/2 0 max(spec)]);
  title('tx spectrum');
  ylabel('Amplitude (dB)')
  xlabel('Frequency (Hz)')

end


% function to write C header file of noise samples so C version gives
% extactly the same results

function write_noise_file(uvnoise_log)

  m = length(uvnoise_log);

  filename = sprintf("../unittest/noise_samples.h");
  f=fopen(filename,"wt");
  fprintf(f,"/* unit variance complex noise samples */\n\n");
  fprintf(f,"/* Generated by write_noise_file() Octave function */\n\n");
  fprintf(f,"COMP noise[]={\n");
  for r=1:m
    if r < m
      fprintf(f, "  {%f,%f},\n", real(uvnoise_log(r)), imag(uvnoise_log(r)));
    else
      fprintf(f, "  {%f,%f}\n};", real(uvnoise_log(r)), imag(uvnoise_log(r)));
    end
  end

  fclose(f);
endfunction
