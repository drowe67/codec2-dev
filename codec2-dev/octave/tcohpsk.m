% tcohpsk.m
% David Rowe Oct 2014
%
% Octave script that tests the C port of the coherent PSK modem.  This
% script loads the output of unittest/tcohpsk.c and compares it to the
% output of the reference versions of the same functions written in
% Octave.
%
% Ideas:
% [ ] EB/No v BER curves changing Np, freq offset etc
%     + can do these pretty fast in C, if we jave the channel models
%     [ ] better interpolation between two signals
%     [ ] feedback to correct out freq offset est
%     [ ] fading channel

graphics_toolkit ("gnuplot");

cohpsk;
fdmdv;
autotest;

rand('state',1); 
randn('state',1);

n = 840;
frames = 35;
framesize = 32;
foff = 0;

EsNodB = 8;
EsNo = 10^(EsNodB/10);
variance = 1/EsNo;

Rs = 50;
Nc = 4;

% --------------------------------------------------------------------------

afdmdv.Fs = 8000;
afdmdv.Nc = Nc-1;
afdmdv.Rs = Rs;
afdmdv.M  = afdmdv.Fs/afdmdv.Rs;
afdmdv.tx_filter_memory = zeros(afdmdv.Nc+1, Nfilter);
afdmdv.Nfilter =  Nfilter;
afdmdv.gt_alpha5_root = gt_alpha5_root;
afdmdv.Fsep = 75;
afdmdv.phase_tx = ones(afdmdv.Nc+1,1);
freq_hz = Fsep*( -Nc/2 - 0.5 + (1:Nc) );
afdmdv.freq_pol = 2*pi*freq_hz/Fs;
afdmdv.freq = exp(j*afdmdv.freq_pol);
afdmdv.Fcentre = 1500;

afdmdv.fbb_rect = exp(j*2*pi*Fcentre/Fs);
afdmdv.fbb_phase_tx = 1;

afdmdv.fbb_rect_rx = exp(j*2*pi*(Fcentre)/Fs);
fbb_phase_rx = 1;

afdmdv.Nrxdec = 31;
afdmdv.rxdec_coeff = fir1(afdmdv.Nrxdec-1, 0.25)';
afdmdv.rxdec_lpf_mem = zeros(1,afdmdv.Nrxdec-1+M);

P = afdmdv.P = 4;
afdmdv.phase_rx = ones(afdmdv.Nc+1,1);
afdmdv.Nsym = 6;
afdmdv.Nfilter = afdmdv.Nsym*afdmdv.M;
afdmdv.rx_fdm_mem = zeros(1,afdmdv.Nfilter + afdmdv.M);
Q = afdmdv.Q = afdmdv.M/4;

afdmdv.Nt = 5;
afdmdv.rx_filter_mem_timing = zeros(afdmdv.Nc+1, afdmdv.Nt*afdmdv.P);
afdmdv.Nfiltertiming = afdmdv.M + afdmdv.Nfilter + afdmdv.M;

afdmdv.rx_filter_memory = zeros(afdmdv.Nc+1, afdmdv.Nfilter);

% ---------------------------------------------------------

load ../build_linux/unittest/tcohpsk_out.txt

acohpsk = standard_init();
acohpsk.framesize        = framesize;
acohpsk.ldpc_code        = 0;
acohpsk.ldpc_code_rate   = 1;
acohpsk.Nc               = Nc;
acohpsk.Rs               = Rs;
acohpsk.Ns               = 4;
acohpsk.Nchip            = 1;
acohpsk.modulation       = 'qpsk';
acohpsk.do_write_pilot_file = 0;
acohpsk = symbol_rate_init(acohpsk);
acohpsk.Ndft = 1024;
acohpsk.f_est = afdmdv.Fcentre;

% -----------------------------------------------------------

rand('state',1); 
tx_bits_coh = round(rand(1,framesize*10));
ptx_bits_coh = 1;

tx_bits_log = [];
tx_symb_log = [];
rx_amp_log = [];
rx_phi_log = [];
ch_symb_log = [];
rx_symb_log = [];
rx_bits_log = [];
tx_bits_prev_log = [];
noise_log = [];
nerr_log = [];
tx_baseband_log = [];
tx_fdm_frame_log = [];
ch_fdm_frame_log = [];
rx_fdm_frame_bb_log = [];

phase = 1;
freq = exp(j*2*pi*foff/acohpsk.Rs);

ch_symb = zeros(acohpsk.Nsymbrowpilot, acohpsk.Nc);

Nerrs = Tbits = 0;

rx_filt_log = [];
rx_fdm_filter_log = [];
rx_baseband_log = [];
rx_fdm_frame_log = [];
f_err_log = [];
f_err_fail = 0;
ct_symb_ff_log = [];

phase_ch = 1;
sync = 0;

prev_tx_bits = [];

% main loop --------------------------------------------------------------------

for i=1:frames
  tx_bits = tx_bits_coh(ptx_bits_coh:ptx_bits_coh+framesize-1);
  ptx_bits_coh += framesize;
  if ptx_bits_coh > length(tx_bits_coh)
    ptx_bits_coh = 1;
  end

  tx_bits_log = [tx_bits_log tx_bits];

  [tx_symb tx_bits prev_tx_sym] = bits_to_qpsk_symbols(acohpsk, tx_bits, [], []);
  tx_symb_log = [tx_symb_log; tx_symb];

  tx_fdm_frame = [];
  for r=1:acohpsk.Nsymbrowpilot
    tx_onesymb = tx_symb(r,:);
    [tx_baseband afdmdv] = tx_filter(afdmdv, tx_onesymb);
    tx_baseband_log = [tx_baseband_log tx_baseband];
    [tx_fdm afdmdv] = fdm_upconvert(afdmdv, tx_baseband);
    tx_fdm_frame = [tx_fdm_frame tx_fdm];
  end
  tx_fdm_frame_log = [tx_fdm_frame_log tx_fdm_frame];

  %
  % Channel --------------------------------------------------------------------
  %

  % [X] calibrated noise for 1% errors
  % [ ] req of x % probability of sync in y ms?
  %      + take foff est, run filter for XXX symbols?
  %      + look for a couple of rows of pilots, take metric
  %      + try to decode frame, say between two pilots?, get corr metric?
  %      + try again
  % [ ] manually test at +/-75 Hz
  % [ ] good BER with freq offset
  % [ ] Integrate offset est into demod
  %    + slow application or block based?
  %    + just need to get a workable first pass for now
  % [ ] module in cohpsk
  % [ ] C port
  % [ ] smaller block size to min ram req
  % [ ] Feeback loop 

  [ch_fdm_frame phase_ch] = freq_shift(tx_fdm_frame, foff, Fs, phase_ch);
  
  % each carrier has power = 2, total power 2Nc, total symbol rate NcRs, noise BW B=Fs
  % Es/No = (C/Rs)/(N/B), N = var = 2NcFs/NcRs(Es/No) = 2Fs/Rs(Es/No)

  variance = 2*Fs/(acohpsk.Rs*EsNo);
  noise = sqrt(variance*0.5)*(randn(1,acohpsk.Nsymbrowpilot*M) + j*randn(1,acohpsk.Nsymbrowpilot*M));
  noise_log = [noise_log noise];

  ch_fdm_frame += noise;

  ch_fdm_frame_log = [ch_fdm_frame_log ch_fdm_frame];

  %
  % Demod ----------------------------------------------------------------------
  %

  % Coarse Freq offset estimation

  next_sync = sync;
  
  [next_sync acohpsk] = coarse_freq_offset_est(acohpsk, afdmdv, ch_fdm_frame, sync, next_sync);

  nin = M;

  % shift entire FDM signal to 0 Hz
 
  afdmdv.fbb_rect_rx = exp(j*2*pi*acohpsk.f_est/afdmdv.Fs);
  rx_fdm_frame_bb = zeros(1, acohpsk.Nsymbrowpilot*M);
  for r=1:acohpsk.Nsymbrowpilot*M
    fbb_phase_rx = fbb_phase_rx*afdmdv.fbb_rect_rx';
    rx_fdm_frame_bb(r) = ch_fdm_frame(r)*fbb_phase_rx;
  end
  mag = abs(fbb_phase_rx);
  fbb_phase_rx /= mag;
  rx_fdm_frame_bb_log = [rx_fdm_frame_bb_log rx_fdm_frame_bb];
  
  % sample rate demod processing

  ch_symb = zeros(acohpsk.Nsymbrowpilot, Nc);
  for r=1:acohpsk.Nsymbrowpilot

    % donwconvert each FDM carrier to Nc separate baseband signals

    [rx_baseband afdmdv] = fdm_downconvert(afdmdv, rx_fdm_frame_bb(1+(r-1)*M:r*M), nin);
    rx_baseband_log = [rx_baseband_log rx_baseband];

    [rx_filt afdmdv] = rx_filter(afdmdv, rx_baseband, nin);

    rx_filt_log = [rx_filt_log rx_filt];

    [rx_onesym rx_timing env afdmdv] = rx_est_timing(afdmdv, rx_filt, nin);     
    ch_symb(r,:) = rx_onesym;
  end
  
  ch_symb_log = [ch_symb_log; ch_symb];

  % coarse timing (frame sync) and initial fine freq est ---------------------------------------------
  
  [next_sync acohpsk] = frame_sync_fine_freq_est(acohpsk, ch_symb, sync, next_sync);
  acohpsk = fine_freq_correct(acohpsk, sync, next_sync);
  ct_symb_ff_log = [ct_symb_ff_log; acohpsk.ct_symb_ff_buf(1:acohpsk.Nsymbrowpilot,:)];

  if (sync == 4) || (next_sync == 4)  
    [rx_symb rx_bits amp_ phi_ EsNo_] = qpsk_symbols_to_bits(acohpsk, acohpsk.ct_symb_ff_buf);
    rx_symb_log = [rx_symb_log; rx_symb];
    rx_amp_log = [rx_amp_log; amp_];
    rx_phi_log = [rx_phi_log; phi_];
    rx_bits_log = [rx_bits_log rx_bits];
    tx_bits_prev_log = [tx_bits_prev_log prev_tx_bits2];

    % BER stats

    error_positions = xor(prev_tx_bits2, rx_bits);
    Nerrs  += sum(error_positions);
    nerr_log = [nerr_log sum(error_positions)];
    Tbits += length(error_positions);
  end

  %printf("f: %d sync: %d next_sync: %d\n", i, sync, next_sync);
  sync = sync_state_machine(sync, next_sync);

  prev_tx_bits2 = prev_tx_bits;
  prev_tx_bits = tx_bits;

end

stem_sig_and_error(1, 111, tx_bits_log_c(1:n), tx_bits_log(1:n) - tx_bits_log_c(1:n), 'tx bits', [1 n -1.5 1.5])
stem_sig_and_error(2, 211, real(tx_symb_log_c(1:n)), real(tx_symb_log(1:n) - tx_symb_log_c(1:n)), 'tx symb re', [1 n -1.5 1.5])
stem_sig_and_error(2, 212, imag(tx_symb_log_c(1:n)), imag(tx_symb_log(1:n) - tx_symb_log_c(1:n)), 'tx symb im', [1 n -1.5 1.5])

stem_sig_and_error(3, 211, real(tx_fdm_log_c(1:n)), real(tx_fdm_frame_log(1:n) - tx_fdm_frame_log_c(1:n)), 'tx fdm frame re', [1 n -10 10])
stem_sig_and_error(3, 212, imag(tx_fdm_log_c(1:n)), imag(tx_fdm_frame_log(1:n) - tx_fdm_frame_log_c(1:n)), 'tx fdm frame im', [1 n -10 10])
stem_sig_and_error(4, 211, real(ch_fdm_frame_log_c(1:n)), real(ch_fdm_frame_log(1:n) - ch_fdm_frame_log_c(1:n)), 'ch fdm frame re', [1 n -10 10])
stem_sig_and_error(4, 212, imag(ch_fdm_frame_log_c(1:n)), imag(ch_fdm_frame_log(1:n) - ch_fdm_frame_log_c(1:n)), 'ch fdm frame im', [1 n -10 10])
stem_sig_and_error(5, 211, real(rx_fdm_frame_bb_log_c(1:n)), real(rx_fdm_frame_bb_log(1:n) - rx_fdm_frame_bb_log_c(1:n)), 'rx fdm frame bb re', [1 n -10 10])
stem_sig_and_error(5, 212, imag(rx_fdm_frame_bb_log_c(1:n)), imag(rx_fdm_frame_bb_log(1:n) - rx_fdm_frame_bb_log_c(1:n)), 'rx fdm frame bb im', [1 n -10 10])

[n m] = size(ch_symb_log);
stem_sig_and_error(6, 211, real(ch_symb_log_c), real(ch_symb_log - ch_symb_log_c), 'ch symb re', [1 n -1.5 1.5])
stem_sig_and_error(6, 212, imag(ch_symb_log_c), imag(ch_symb_log - ch_symb_log_c), 'ch symb im', [1 n -1.5 1.5])
stem_sig_and_error(7, 211, real(ct_symb_ff_log_c), real(ct_symb_ff_log - ct_symb_ff_log_c), 'ct symb ff re', [1 n -1.5 1.5])
stem_sig_and_error(7, 212, imag(ct_symb_ff_log_c), imag(ct_symb_ff_log - ct_symb_ff_log_c), 'ct symb ff im', [1 n -1.5 1.5])

stem_sig_and_error(8, 211, rx_amp_log_c, rx_amp_log - rx_amp_log_c, 'Amp Est', [1 n -1.5 1.5])
stem_sig_and_error(8, 212, rx_phi_log_c, rx_phi_log - rx_phi_log_c, 'Phase Est', [1 n -4 4])
stem_sig_and_error(9, 211, real(rx_symb_log_c), real(rx_symb_log - rx_symb_log_c), 'rx symb re', [1 n -1.5 1.5])
stem_sig_and_error(9, 212, imag(rx_symb_log_c), imag(rx_symb_log - rx_symb_log_c), 'rx symb im', [1 n -1.5 1.5])
stem_sig_and_error(10, 111, rx_bits_log_c, rx_bits_log - rx_bits_log_c, 'rx bits', [1 n -1.5 1.5])

check(tx_bits_log, tx_bits_log_c, 'tx_bits');
check(tx_symb_log, tx_symb_log_c, 'tx_symb');
check(tx_fdm_frame_log, tx_fdm_frame_log_c, 'tx_fdm_frame');
check(ch_fdm_frame_log, ch_fdm_frame_log_c, 'ch_fdm_frame');
check(rx_fdm_frame_bb_log, rx_fdm_frame_bb_log_c, 'rx_fdm_frame_bb', 0.01);
check(ch_symb_log, ch_symb_log_c, 'ch_symb',0.01);
check(ct_symb_ff_log, ct_symb_ff_log_c, 'ct_symb_ff',0.01);
check(rx_amp_log, rx_amp_log_c, 'rx_amp_log',0.01);
check(rx_phi_log, rx_phi_log_c, 'rx_phi_log',0.01);
check(rx_symb_log, rx_symb_log_c, 'rx_symb',0.01);
check(rx_bits_log, rx_bits_log_c, 'rx_bits');

% Determine bit error rate

sz = length(tx_bits_log_c);
Nerrs_c = sum(xor(tx_bits_prev_log, rx_bits_log_c));
Tbits_c = length(tx_bits_prev_log);
ber_c = Nerrs_c/Tbits_c;
ber = Nerrs/Tbits;
printf("EsNodB: %4.1f ber..: %3.2f Nerrs..: %d Tbits..: %d\n", EsNodB, ber, Nerrs, Tbits);
printf("EsNodB: %4.1f ber_c: %3.2f Nerrs_c: %d Tbits_c: %d\n", EsNodB, ber_c, Nerrs_c, Tbits_c);

% C header file of noise samples so C version gives extacly the same results

function write_noise_file(noise_log)

  m = length(noise_log);

  filename = sprintf("../unittest/noise_samples.h");
  f=fopen(filename,"wt");
  fprintf(f,"/* Generated by write_noise_file() Octave function */\n\n");
  fprintf(f,"COMP noise[]={\n");
  for r=1:m
    if r < m
      fprintf(f, "  {%f,%f},\n", real(noise_log(r)), imag(noise_log(r)));
    else
      fprintf(f, "  {%f,%f}\n};", real(noise_log(r)), imag(noise_log(r)));
    end
  end

  fclose(f);
endfunction
