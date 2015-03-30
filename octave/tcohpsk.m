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

n = 2000;
frames = 35*4;
framesize = 32;
foff = 0;

EsNodB = 8;
EsNo = 10^(EsNodB/10);
variance = 1/EsNo;

Rs = 50
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
noise_log = [];
nerr_log = [];
tx_baseband_log = [];
tx_fdm_log = [];

phase = 1;
freq = exp(j*2*pi*foff/acohpsk.Rs);

ch_symb = zeros(acohpsk.Nsymbrowpilot, acohpsk.Nc);

Nerrs = Tbits = 0;

rx_filt_log = [];
rx_fdm_filter_log = [];
rx_baseband_log = [];
rx_fdm_log = [];
f_err_log = [];
f_err_fail = 0;

fbb_phase_ch = 1;
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
  tx_fdm_log = [tx_fdm_log tx_fdm_frame];

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

  afdmdv.fbb_rect_ch = exp(j*2*pi*foff/Fs);
  ch_fdm_frame = zeros(1, acohpsk.Nsymbrowpilot*M);
  for r=1:acohpsk.Nsymbrowpilot*M
    fbb_phase_ch = fbb_phase_ch*afdmdv.fbb_rect_ch;
    ch_fdm_frame(r) = tx_fdm_frame(r)*fbb_phase_ch;
  end
  mag = abs(fbb_phase_ch);
  fbb_phase_ch /= mag;
  
  % each carrier has power = 2, total power 2Nc, total symbol rate NcRs, noise BW B=Fs
  % Es/No = (C/Rs)/(N/B), N = var = 2NcFs/NcRs(Es/No) = 2Fs/Rs(Es/No)

  variance = 2*Fs/(acohpsk.Rs*EsNo);
  noise = sqrt(variance*0.5)*(randn(1,acohpsk.Nsymbrowpilot*M) + j*randn(1,acohpsk.Nsymbrowpilot*M));
  noise_log = [noise_log; noise];

  ch_fdm_frame += noise;

  %
  % Demod ----------------------------------------------------------------------
  %

  % Coarse Freq offset estimation

  next_sync = sync;
  
  [next_sync acohpsk] = coarse_freq_offset_est(acohpsk, afdmdv, ch_fdm_frame, sync, next_sync);

  nin = M;

  % shift frame down to complex baseband
 
  afdmdv.fbb_rect_rx = exp(j*2*pi*acohpsk.f_est/afdmdv.Fs);
  rx_fdm_frame_bb = zeros(1, acohpsk.Nsymbrowpilot*M);
  for r=1:acohpsk.Nsymbrowpilot*M
    fbb_phase_rx = fbb_phase_rx*afdmdv.fbb_rect_rx';
    rx_fdm_frame_bb(r) = ch_fdm_frame(r)*fbb_phase_rx;
  end
  mag = abs(fbb_phase_rx);
  fbb_phase_rx /= mag;
  rx_fdm_log = [rx_fdm_log rx_fdm_frame_bb];
  
  % sample rate demod processing

  ch_symb = zeros(acohpsk.Nsymbrowpilot, Nc);
  for r=1:acohpsk.Nsymbrowpilot

    [rx_baseband afdmdv] = fdm_downconvert(afdmdv, rx_fdm_frame_bb(1+(r-1)*M:r*M), nin);
    rx_baseband_log = [rx_baseband_log rx_baseband];
    [rx_filt afdmdv] = rx_filter(afdmdv, rx_baseband, nin);

    rx_filt_log = [rx_filt_log rx_filt];

    [rx_onesym rx_timing env afdmdv] = rx_est_timing(afdmdv, rx_filt, nin);     
    ch_symb(r,:) = rx_onesym;
  end
  
  ch_symb_log = [ch_symb_log; ch_symb];

  % coarse timing (frame sync) and initial fine freq est ---------------------------------------------
  
  [next_sync acohpsk] = frame_sync_fine_timing_est(acohpsk, ch_symb, sync, next_sync);
  %acohpsk.ff_rect = exp(j*2*pi*(-1.0)/Rs);

  if (i==1000)
    xx
  end

  acohpsk = fine_freq_correct(acohpsk, sync, next_sync);

  if (sync == 4) || (next_sync == 4)  
    [rx_symb rx_bits amp_ phi_ EsNo_] = qpsk_symbols_to_bits(acohpsk, acohpsk.ct_symb_ff_buf);
    rx_symb_log = [rx_symb_log; rx_symb];
    rx_amp_log = [rx_amp_log; amp_];
    rx_phi_log = [rx_phi_log; phi_];
    rx_bits_log = [rx_bits_log rx_bits];

    % BER stats

    if i > 2
      error_positions = xor(prev_tx_bits2, rx_bits);
      Nerrs  += sum(error_positions);
      nerr_log = [nerr_log sum(error_positions)];
     Tbits += length(error_positions);
    end 
  end

  sync = sync_state_machine(sync, next_sync);

  prev_tx_bits2 = prev_tx_bits;
  prev_tx_bits = tx_bits;

end

stem_sig_and_error(1, 111, tx_bits_log_c(1:n), tx_bits_log(1:n) - tx_bits_log_c(1:n), 'tx bits', [1 n -1.5 1.5])
stem_sig_and_error(2, 211, real(tx_symb_log_c(1:n)), real(tx_symb_log(1:n) - tx_symb_log_c(1:n)), 'tx symb re', [1 n -1.5 1.5])
stem_sig_and_error(2, 212, imag(tx_symb_log_c(1:n)), imag(tx_symb_log(1:n) - tx_symb_log_c(1:n)), 'tx symb im', [1 n -1.5 1.5])

stem_sig_and_error(3, 211, real(tx_fdm_log_c(1:n)), real(tx_fdm_log(1:n) - tx_fdm_log_c(1:n)), 'tx fdm re', [1 n -10 10])
stem_sig_and_error(3, 212, imag(tx_fdm_log_c(1:n)), imag(tx_fdm_log(1:n) - tx_fdm_log_c(1:n)), 'tx fdm im', [1 n -10 10])

stem_sig_and_error(4, 211, real(ch_symb_log_c(1:n)), real(ch_symb_log(1:n) - ch_symb_log_c(1:n)), 'ch symb re', [1 n -2 2])
stem_sig_and_error(4, 212, imag(ch_symb_log_c(1:n)), imag(ch_symb_log(1:n) - ch_symb_log_c(1:n)), 'ch symb im', [1 n -2 2])
stem_sig_and_error(5, 211, rx_amp_log_c(1:n), rx_amp_log(1:n) - rx_amp_log_c(1:n), 'Amp Est', [1 n -1.5 1.5])
stem_sig_and_error(5, 212, rx_phi_log_c(1:n), rx_phi_log(1:n) - rx_phi_log_c(1:n), 'Phase Est', [1 n -4 4])
stem_sig_and_error(6, 211, real(rx_symb_log_c(1:n)), real(rx_symb_log(1:n) - rx_symb_log_c(1:n)), 'rx symb re', [1 n -1.5 1.5])
stem_sig_and_error(6, 212, imag(rx_symb_log_c(1:n)), imag(rx_symb_log(1:n) - rx_symb_log_c(1:n)), 'rx symb im', [1 n -1.5 1.5])
stem_sig_and_error(7, 111, rx_bits_log_c(1:n), rx_bits_log(1:n) - rx_bits_log_c(1:n), 'rx bits', [1 n -1.5 1.5])

check(tx_bits_log, tx_bits_log_c, 'tx_bits');
check(tx_symb_log, tx_symb_log_c, 'tx_symb');
check(tx_fdm_log, tx_fdm_log_c, 'tx_fdm');
check(rx_fdm_log, rx_fdm_log_c, 'rx_fdm');
if 0
check(rx_baseband_log, rx_baseband_log_c, 'rx_baseband',0.01);
check(rx_filt_log, rx_filt_log_c, 'rx_filt');
check(ch_symb_log, ch_symb_log_c, 'ch_symb',0.01);
check(rx_amp_log, rx_amp_log_c, 'rx_amp_log',0.01);
check(rx_phi_log, rx_phi_log_c, 'rx_phi_log');
check(rx_symb_log, rx_symb_log_c, 'rx_symb',0.01);
check(rx_bits_log, rx_bits_log_c, 'rx_bits');
end

% Determine bit error rate

sz = length(tx_bits_log_c);
Nerrs_c = sum(xor(tx_bits_log_c(framesize+1:sz-2*framesize), rx_bits_log_c(3*framesize+1:sz)));
Tbits_c = sz - 2*framesize;
ber_c = Nerrs_c/Tbits_c;
ber = Nerrs/Tbits;
printf("EsNodB: %4.1f ber..: %3.2f Nerrs..: %d Tbits..: %d\n", EsNodB, ber, Nerrs, Tbits);
printf("EsNodB: %4.1f ber_c: %3.2f Nerrs_c: %d Tbits_c: %d\n", EsNodB, ber_c, Nerrs_c, Tbits_c);
printf("f_err std: %f  fails: %d\n", std(f_err_log), f_err_fail);

% C header file of noise samples so C version gives extacly the same results

function write_noise_file(noise_log)

  [m n] = size(noise_log);

  filename = sprintf("../unittest/noise_samples.h");
  f=fopen(filename,"wt");
  fprintf(f,"/* Generated by write_noise_file() Octave function */\n\n");
  fprintf(f,"COMP noise[][PILOTS_NC]={\n");
  for r=1:m
    fprintf(f, "  {");
    for c=1:n-1
      fprintf(f, "  {%f,%f},", real(noise_log(r, c)), imag(noise_log(r, c)));
    end
    if r < m
      fprintf(f, "  {%f,%f}},\n", real(noise_log(r, n)), imag(noise_log(r, n)));
    else
      fprintf(f, "  {%f,%f}}\n};", real(noise_log(r, n)), imag(noise_log(r, n)));
    end
  end

  fclose(f);
endfunction
