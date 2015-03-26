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
frames = 35;
framesize = 160
foff = 50;

EsNodB = 18;
EsNo = 10^(EsNodB/10);
variance = 1/EsNo;

load ../build_linux/unittest/tcohpsk_out.txt

sim_in = standard_init();
sim_in.framesize        = 160;
sim_in.ldpc_code        = 0;
sim_in.ldpc_code_rate   = 1;
sim_in.Nc = Nc          = 4;
sim_in.Rs               = 50;
sim_in.Ns               = 4;
sim_in.Np               = 2;
sim_in.Nchip            = 1;
sim_in.modulation       = 'qpsk';
sim_in.do_write_pilot_file = 0;
sim_in = symbol_rate_init(sim_in);

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
freq = exp(j*2*pi*foff/sim_in.Rs);

ch_symb = zeros(sim_in.Nsymbrowpilot, sim_in.Nc);

Nerrs = Tbits = 0;

fdmdv.Fs = 8000;
fdmdv.Nc = Nc-1;
M = fdmdv.M = Fs/Rs;
fdmdv.tx_filter_memory = zeros(fdmdv.Nc+1, Nfilter);
fdmdv.Nfilter =  Nfilter;
fdmdv.gt_alpha5_root = gt_alpha5_root;
fdmdv.Fsep = 75;
fdmdv.phase_tx = ones(fdmdv.Nc+1,1);
freq_hz = Fsep*( -Nc/2 - 0.5 + (1:Nc) );
fdmdv.freq_pol = 2*pi*freq_hz/Fs;
fdmdv.freq = exp(j*fdmdv.freq_pol);
Fcentre = 1500;

fdmdv.fbb_rect = exp(j*2*pi*Fcentre/Fs);
fdmdv.fbb_phase_tx = 1;

fdmdv.fbb_rect_rx = exp(j*2*pi*(Fcentre)/Fs);
fbb_phase_rx = 1;

fdmdv.Nrxdec = 31;
fdmdv.rxdec_coeff = fir1(fdmdv.Nrxdec-1, 0.25)';
fdmdv.rxdec_lpf_mem = zeros(1,fdmdv.Nrxdec-1+M);

P = fdmdv.P = 4;
fdmdv.phase_rx = ones(fdmdv.Nc+1,1);
fdmdv.Nsym = 6;
fdmdv.Nfilter = fdmdv.Nsym*fdmdv.M;
fdmdv.rx_fdm_mem = zeros(1,fdmdv.Nfilter + fdmdv.M);
Q = fdmdv.Q = fdmdv.M/4;

fdmdv.Nt = 5;
fdmdv.rx_filter_mem_timing = zeros(fdmdv.Nc+1, fdmdv.Nt*fdmdv.P);
fdmdv.Nfiltertiming = fdmdv.M + fdmdv.Nfilter + fdmdv.M;

fdmdv.rx_filter_memory = zeros(fdmdv.Nc+1, fdmdv.Nfilter);

rx_filt_log = [];
rx_fdm_filter_log = [];
rx_baseband_log = [];
rx_fdm_log = [];
f_err_log = [];
f_err_fail = 0;

fbb_phase_ch = 1;
sync = 0;

% set up pilot signals for sync ------------------------------------------------
 
% frame of just pilots for coarse sync

tx_bits = zeros(1,framesize);
[tx_symb_pilot tx_bits prev_tx_sym] = bits_to_qpsk_symbols(sim_in, tx_bits, [], []);
for r=1:(sim_in.Ns+1):sim_in.Nsymbrowpilot
  tx_symb_pilot(r+1:r+sim_in.Ns,:) = zeros(sim_in.Ns, sim_in.Nc);
end

% filtered frame of just pilots for freq offset and coarse timing

tx_pilot_fdm_frame = [];
fdmdvp = fdmdv;
for r=1:sim_in.Nsymbrowpilot
  tx_onesymb = tx_symb_pilot(r,:);
  [tx_baseband fdmdvp] = tx_filter(fdmdvp, tx_onesymb);
  tx_baseband_log = [tx_baseband_log tx_baseband];
  [tx_fdm fdmdvp] = fdm_upconvert(fdmdvp, tx_baseband);
  tx_pilot_fdm_frame = [tx_pilot_fdm_frame tx_fdm];
end

% ------------------------------------------------------------------------------

ct_symb_buf = zeros(2*sim_in.Nsymbrowpilot, sim_in.Nc);

prev_tx_bits = [];

% main loop --------------------------------------------------------------------

for i=1:frames
  tx_bits = tx_bits_coh(ptx_bits_coh:ptx_bits_coh+framesize-1);
  ptx_bits_coh += framesize;
  if ptx_bits_coh > length(tx_bits_coh)
    ptx_bits_coh = 1;
  end

  tx_bits_log = [tx_bits_log tx_bits];

  [tx_symb tx_bits prev_tx_sym] = bits_to_qpsk_symbols(sim_in, tx_bits, [], []);
  tx_symb_log = [tx_symb_log; tx_symb];

  tx_fdm_frame = [];
  for r=1:sim_in.Nsymbrowpilot
    tx_onesymb = tx_symb(r,:);
    [tx_baseband fdmdv] = tx_filter(fdmdv, tx_onesymb);
    tx_baseband_log = [tx_baseband_log tx_baseband];
    [tx_fdm fdmdv] = fdm_upconvert(fdmdv, tx_baseband);
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

  fdmdv.fbb_rect_ch = exp(j*2*pi*foff/Fs);
  ch_fdm_frame = zeros(1, sim_in.Nsymbrowpilot*M);
  for r=1:sim_in.Nsymbrowpilot*M
    fbb_phase_ch = fbb_phase_ch*fdmdv.fbb_rect_ch;
    ch_fdm_frame(r) = tx_fdm_frame(r)*fbb_phase_ch;
  end
  mag = abs(fbb_phase_ch);
  fbb_phase_ch /= mag;
  
  % each carrier has power = 2, total power 2Nc, total symbol rate NcRs, noise BW B=Fs
  % Es/No = (C/Rs)/(N/B), N = var = 2NcFs/NcRs(Es/No) = 2Fs/Rs(Es/No)

  variance = 2*Fs/(sim_in.Rs*EsNo);
  noise = sqrt(variance*0.5)*(randn(1,sim_in.Nsymbrowpilot*M) + j*randn(1,sim_in.Nsymbrowpilot*M));
  noise_log = [noise_log; noise];

  ch_fdm_frame += noise;

  %
  % Demod ----------------------------------------------------------------------
  %

  % Coarse Freq offset estimation

  if sync == 0
    f_start = Fcentre - ((Nc/2)+2)*Fsep; f_stop = Fcentre + ((Nc/2)+2)*Fsep;
    T = abs(fft(ch_fdm_frame(1:8*M).* hanning(8*M)',Fs)).^2;
    T  = T(f_start:f_stop);
    x = f_start:f_stop;
    f_est = x*T'/sum(T);
    f_off_est = f_est - Fcentre;
    f_err = f_off_est - foff;
    if abs(f_err) > 5
      f_err_fail++;
    end
    f_err_log = [f_err_log f_err];
    printf("coarse freq est: %f offset: %f  f_err: %f\n", f_est, f_off_est, f_err);
    fdmdv.fbb_rect_rx = exp(j*2*pi*(f_est)/Fs);
    sync = 1;
  end

  nin = M;

  % shift frame down to complex baseband

  rx_fdm_frame_bb = zeros(1, sim_in.Nsymbrowpilot*M);
  for r=1:sim_in.Nsymbrowpilot*M
    fbb_phase_rx = fbb_phase_rx*fdmdv.fbb_rect_rx';
    rx_fdm_frame_bb(r) = ch_fdm_frame(r)*fbb_phase_rx;
  end
  mag = abs(fbb_phase_rx);
  fbb_phase_rx /= mag;
  rx_fdm_log = [rx_fdm_log rx_fdm_frame_bb];
  
  ch_symb = zeros(sim_in.Nsymbrowpilot, Nc);
  for r=1:sim_in.Nsymbrowpilot

    [rx_baseband fdmdv] = fdm_downconvert(fdmdv, rx_fdm_frame_bb(1+(r-1)*M:r*M), nin);
    rx_baseband_log = [rx_baseband_log rx_baseband];
    [rx_filt fdmdv] = rx_filter(fdmdv, rx_baseband, nin);

    rx_filt_log = [rx_filt_log rx_filt];

    [rx_onesym rx_timing env fdmdv] = rx_est_timing(fdmdv, rx_filt, nin);     
    ch_symb(r,:) = rx_onesym;
  end
  
  ch_symb_log = [ch_symb_log; ch_symb];

  % coarse timing (frame sync) and initial fine freq est ---------------------------------------------

  ct_symb_buf(1:sim_in.Nsymbrowpilot,:) = ct_symb_buf(sim_in.Nsymbrowpilot+1:2*sim_in.Nsymbrowpilot,:);
  ct_symb_buf(sim_in.Nsymbrowpilot+1:2*sim_in.Nsymbrowpilot,:) = ch_symb;

  sampling_points = [1 (2:sim_in.Ns+1:1+sim_in.Npilotsframe*sim_in.Ns+1)];
  pilot2 = [ sim_in.pilot(1,:); sim_in.pilot];

  if sync == 1
    max_corr = 0;
    for f_fine=-15:1:15
      f_fine_rect = exp(-j*f_fine*2*pi*sampling_points/Rs)';
      for t=0:sim_in.Nsymbrowpilot-1
        corr = 0; mag = 0;
        for c=1:Nc
          f_corr_vec = f_fine_rect .* ct_symb_buf(t+sampling_points,c);
          for p=1:sim_in.Npilotsframe+1
            corr += pilot2(p,c) * f_corr_vec(p);
            mag  += abs(f_corr_vec(p));
          end
        end
        %printf("  f: %f  t: %d corr: %f %f\n", f_fine, t, real(corr), imag(corr));
        if corr > max_corr
          max_corr = corr;
          max_mag = mag;
          ct = t;
          f_fine_est = f_fine;
        end
      end
    end

    printf("  fine freq f: %f max_corr: %f max_mag: %f ct: %d\n", f_fine_est, abs(max_corr), max_mag, ct);
    if max_corr/max_mag > 0.8
      sync = 2;
    end
  end
 
  if (i==50)
    figure(8)
    f_fine_rect = exp(-j*f_fine_est*2*pi*sampling_points/Rs)';
    plot(f_fine_rect,'+');
    hold on;
    plot(ct_symb_buf(ct+sampling_points,1),'b+');
    hold off; 
    xx
  end

  [rx_symb rx_bits rx_symb_linear amp_linear amp_ phi_ EsNo_ prev_sym_rx sim_in] = qpsk_symbols_to_bits(sim_in, ct_symb_buf(ct+1:ct+sim_in.Nsymbrowpilot,:), []);
  rx_symb_log = [rx_symb_log; rx_symb];
  rx_amp_log = [rx_amp_log; amp_];
  rx_phi_log = [rx_phi_log; phi_];
  rx_bits_log = [rx_bits_log rx_bits];

  % BER stats

  if i > 3
    error_positions = xor(prev_tx_bits2, rx_bits);
    Nerrs  += sum(error_positions);
    nerr_log = [nerr_log sum(error_positions)];
    Tbits += length(error_positions);
  end
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
check(rx_baseband_log, rx_baseband_log_c, 'rx_baseband',0.01);
check(rx_filt_log, rx_filt_log_c, 'rx_filt');
check(ch_symb_log, ch_symb_log_c, 'ch_symb',0.01);
check(rx_amp_log, rx_amp_log_c, 'rx_amp_log',0.01);
check(rx_phi_log, rx_phi_log_c, 'rx_phi_log');
check(rx_symb_log, rx_symb_log_c, 'rx_symb',0.01);
check(rx_bits_log, rx_bits_log_c, 'rx_bits');

% Determine bit error rate

sz = length(tx_bits_log_c);
Nerrs_c = sum(xor(tx_bits_log_c(framesize+1:sz-2*framesize), rx_bits_log_c(3*framesize+1:sz)));
Tbits_c = sz - 2*framesize;
ber_c = Nerrs_c/Tbits_c;
ber = Nerrs/Tbits;
printf("EsNodB: %4.1f ber..: %3.2f Nerrs..: %d Tbits..: %d\n", EsNodB, ber, Nerrs, Tbits);
printf("EsNodB: %4.1f ber_c: %3.2f Nerrs_c: %d Tbits_c: %d\n", EsNodB, ber_c, Nerrs_c, Tbits_c);
printf("f_err std: %f  fails: %d\n", std(f_err_log), f_err_fail);
figure(8)
hist(f_err_log)

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
