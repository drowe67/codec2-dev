% test_foff.m
% David Rowe April 2015
%
% Octave script for testing the cohpsk freq offset estimator

graphics_toolkit ("gnuplot");
more off;

cohpsk;
fdmdv;
autotest;

rand('state',1); 
randn('state',1);

 
% Core function for testing frequency offset estimator.  Performs one test

function sim_out = freq_off_est_test(sim_in)
  global Nfilter;
  global M;

  Rs = 100;
  Nc = 4;
  Nd = 2;
  framesize = 32;
  Fs = 8000;
  Fcentre = 1500;
  
  Nsw = 3; % numbers of sync window frames we process over.  Set based
           % on Nsym to flush filter memory by final frame in windw

  afdmdv.Nsym = 5;
  afdmdv.Nt = 5;

  afdmdv.Fs = 8000;
  afdmdv.Nc = Nd*Nc-1;
  afdmdv.Rs = Rs;
  M = afdmdv.M  = afdmdv.Fs/afdmdv.Rs;
  afdmdv.Nfilter = afdmdv.Nsym*M;
  afdmdv.tx_filter_memory = zeros(afdmdv.Nc+1, afdmdv.Nfilter);
  excess_bw = 0.5;
  afdmdv.gt_alpha5_root = gen_rn_coeffs(excess_bw, 1/Fs, Rs, afdmdv.Nsym, afdmdv.M);
  afdmdv.Fsep = afdmdv.Rs*(1+excess_bw);
  afdmdv.phase_tx = ones(afdmdv.Nc+1,1);
  freq_hz = afdmdv.Fsep*( -Nc*Nd/2 - 0.5 + (1:Nc*Nd).^1.1 );
  afdmdv.freq_pol = 2*pi*freq_hz/Fs;
  afdmdv.freq = exp(j*afdmdv.freq_pol);
  afdmdv.Fcentre = 1500;

  afdmdv.fbb_rect = exp(j*2*pi*Fcentre/Fs);
  afdmdv.fbb_phase_tx = 1;
  afdmdv.fbb_phase_rx = 1;
  %afdmdv.phase_rx = ones(afdmdv.Nc+1,1);
  afdmdv.phase_rx = 1 - 2*(mod(1:Nc*Nd,2));
  nin = M;

  P = afdmdv.P = 4;
  afdmdv.Nfilter = afdmdv.Nsym*afdmdv.M;
  afdmdv.rx_filter_mem_timing = zeros(afdmdv.Nc+1, afdmdv.Nt*afdmdv.P);
  afdmdv.Nfiltertiming = afdmdv.M + afdmdv.Nfilter + afdmdv.M;
  afdmdv.rx_filter_memory = zeros(afdmdv.Nc+1, afdmdv.Nfilter);

  acohpsk = standard_init();
  acohpsk.framesize        = framesize;
  acohpsk.ldpc_code        = 0;
  acohpsk.ldpc_code_rate   = 1;
  acohpsk.Nc               = Nc;
  acohpsk.Rs               = Rs;
  acohpsk.Ns               = 4;
  acohpsk.Nd               = Nd;
  acohpsk.modulation       = 'qpsk';
  acohpsk.do_write_pilot_file = 0;
  acohpsk = symbol_rate_init(acohpsk);
  acohpsk.Ncm  = 10*acohpsk.Nsymbrowpilot*M;
  acohpsk.coarse_mem  = zeros(1,acohpsk.Ncm);
  acohpsk.Ndft = 2^(ceil(log2(acohpsk.Ncm)));
 
  ch_fdm_frame_buf = zeros(1, Nsw*acohpsk.Nsymbrowpilot*M);

  frames    = sim_in.frames;
  EsNodB    = sim_in.EsNodB;
  foff      = sim_in.foff;
  dfoff     = sim_in.dfoff;
  fading_en = sim_in.fading_en;

  EsNo = 10^(EsNodB/10);
  hf_delay_ms = 2;
  phase_ch = 1;

  rand('state',1); 
  tx_bits_coh = round(rand(1,framesize*10));
  ptx_bits_coh = 1;
  [spread spread_2ms hf_gain] = init_hf_model(Fs, frames*acohpsk.Nsymbrowpilot*afdmdv.M);

  hf_n = 1;
  nhfdelay = floor(hf_delay_ms*Fs/1000);
  ch_fdm_delay = zeros(1, acohpsk.Nsymbrowpilot*M + nhfdelay);
  
  sync = 0; next_sync = 1;
  sync_start = 1;
  freq_offset_log = [];
  sync_time_log = [];
  ch_fdm_frame_log = [];
  ch_symb_log = [];
  tx_fdm_frame_log = [];
  ratio_log = [];

  for f=1:frames

    acohpsk.frame = f;

    %
    % Mod --------------------------------------------------------------------
    %

    tx_bits = tx_bits_coh(ptx_bits_coh:ptx_bits_coh+framesize-1);
    ptx_bits_coh += framesize;
    if ptx_bits_coh > length(tx_bits_coh)
      ptx_bits_coh = 1;
    end 

    [tx_symb tx_bits] = bits_to_qpsk_symbols(acohpsk, tx_bits, [], []);

    tx_fdm_frame = [];
    for r=1:acohpsk.Nsymbrowpilot
      tx_onesymb = tx_symb(r,:);
      [tx_baseband afdmdv] = tx_filter(afdmdv, tx_onesymb);
      [tx_fdm afdmdv] = fdm_upconvert(afdmdv, tx_baseband);
      tx_fdm_frame = [tx_fdm_frame tx_fdm];
    end
    tx_fdm_frame_log = [tx_fdm_frame_log tx_fdm_frame];

    %
    % Channel --------------------------------------------------------------------
    %

    ch_fdm_frame = zeros(1,acohpsk.Nsymbrowpilot*M);
    for i=1:acohpsk.Nsymbrowpilot*M
      foff_rect = exp(j*2*pi*foff/Fs);
      foff += dfoff;
      phase_ch *= foff_rect;
      ch_fdm_frame(i) = tx_fdm_frame(i) * phase_ch;
    end
    phase_ch /= abs(phase_ch);

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
    noise = sqrt(variance)*uvnoise;

    ch_fdm_frame += noise;
    ch_fdm_frame_log = [ch_fdm_frame_log ch_fdm_frame];

    %
    % Try to achieve sync --------------------------------------------------------------------
    %

    % store Nsw frames so we can rewind if we get a good candidate

    ch_fdm_frame_buf(1:(Nsw-1)*acohpsk.Nsymbrowpilot*M) = ch_fdm_frame_buf(acohpsk.Nsymbrowpilot*M+1:Nsw*acohpsk.Nsymbrowpilot*M);
    ch_fdm_frame_buf((Nsw-1)*acohpsk.Nsymbrowpilot*M+1:Nsw*acohpsk.Nsymbrowpilot*M) = ch_fdm_frame;

    next_sync = sync;
    sync = 0;

    if sync == 0

      % we can test +/- 20Hz, so we break this up into 3 tests to cover +/- 60Hz

      max_ratio = 0;
      for acohpsk.f_est = Fcentre-40:40:Fcentre+40
      %for acohpsk.f_est = Fcentre
        
        printf("  [%d] acohpsk.f_est: %f +/- 20\n", f, acohpsk.f_est);
        [ch_symb rx_timing rx_filt rx_baseband afdmdv acohpsk.f_est ] = rate_Fs_rx_processing(afdmdv, ch_fdm_frame_buf, acohpsk.f_est, Nsw*acohpsk.Nsymbrowpilot, nin, 0);

        % coarse timing (frame sync) and initial fine freq est ---------------------------------------------
  
        for i=1:Nsw-1
          acohpsk.ct_symb_buf = update_ct_symb_buf(acohpsk.ct_symb_buf, ch_symb((i-1)*acohpsk.Nsymbrowpilot+1:i*acohpsk.Nsymbrowpilot,:), acohpsk.Nct_sym_buf, acohpsk.Nsymbrowpilot);
        end
        [anext_sync acohpsk] = frame_sync_fine_freq_est(acohpsk, ch_symb((Nsw-1)*acohpsk.Nsymbrowpilot+1:Nsw*acohpsk.Nsymbrowpilot,:), sync, next_sync);

        if anext_sync == 1
          printf("  [%d] acohpsk.ratio: %f\n", f, acohpsk.ratio);
          if acohpsk.ratio > max_ratio
            max_ratio   = acohpsk.ratio;
            f_est       = acohpsk.f_est - acohpsk.f_fine_est;
            next_sync   = anext_sync;
          end
        end
      end
    end

    % we've found a sync candidate

    if (next_sync == 1)

       % rewind and re-process last few frames with f_est

       acohpsk.f_est = f_est;
       printf("  [%d] trying sync and f_est: %f\n", f, acohpsk.f_est);
       [ch_symb rx_timing rx_filt rx_baseband afdmdv acohpsk.f_est] = rate_Fs_rx_processing(afdmdv, ch_fdm_frame_buf, acohpsk.f_est, Nsw*acohpsk.Nsymbrowpilot, nin, 0);
       ch_symb_log = ch_symb;
       for i=1:Nsw-1
         acohpsk.ct_symb_buf = update_ct_symb_buf(acohpsk.ct_symb_buf, ch_symb((i-1)*acohpsk.Nsymbrowpilot+1:i*acohpsk.Nsymbrowpilot,:), acohpsk.Nct_sym_buf, acohpsk.Nsymbrowpilot);
       end
       [next_sync acohpsk] = frame_sync_fine_freq_est(acohpsk, ch_symb((Nsw-1)*acohpsk.Nsymbrowpilot+1:Nsw*acohpsk.Nsymbrowpilot,:), sync, next_sync);
       if abs(acohpsk.f_fine_est) > 2
         printf("  [%d] Hmm %f is a bit big so back to coarse est ...\n", f, acohpsk.f_fine_est);
         next_sync = 0;
       end

       % candidate checks out so log stats

       if (next_sync == 1)
         printf("  [%d] in sync!\n", f);
         freq_offset_log = [freq_offset_log Fcentre+foff-acohpsk.f_est];
         sync_time_log = [sync_time_log f-sync_start];
         ratio_log = [ratio_log max_ratio];
         next_sync = 0; sync_start = f;
       end
    end

    %printf("f: %d sync: %d next_sync: %d\n", f, sync, next_sync);
    [sync acohpsk] = sync_state_machine(acohpsk, sync, next_sync);

  end

  % ftx=fopen("coarse_tx.raw","wb"); fwrite(ftx, 1000*ch_fdm_frame_log, "short"); fclose(ftx);

  sim_out.freq_offset_log = freq_offset_log;
  sim_out.sync_time_log = sync_time_log;
  sim_out.ch_fdm_frame_log = ch_fdm_frame_log;
  sim_out.ch_symb_log = ch_symb_log;
  sim_out.tx_fdm_frame_log = tx_fdm_frame_log;
  sim_out.ratio_log = ratio_log;
endfunction


function freq_off_est_test_single
  sim_in.frames    = 5;
  sim_in.EsNodB    = 120;
  sim_in.foff      = 0;
  sim_in.dfoff     = 0;
  sim_in.fading_en = 0;

  sim_out = freq_off_est_test(sim_in);

  figure(1)
  subplot(211)
  plot(sim_out.freq_offset_log)
  ylabel('f est error')
  xlabel('time')

  subplot(212)
  if length(sim_out.freq_offset_log) > 0
    hist(sim_out.freq_offset_log)
    xlabel('f est error')
    ylabel('histogram');
  end

  figure(2)
  subplot(211)
  plot(sim_out.sync_time_log)
  ylabel('time to sync (frames)')
  subplot(212)
  if length(sim_out.sync_time_log) > 0
    hist(sim_out.sync_time_log)
    xlabel('time to sync (frames)')
    ylabel('histogram')
  end

  figure(3)
  subplot(211)
  plot(real(sim_out.tx_fdm_frame_log(1:2*960)))
  grid;
  subplot(212)
  plot(real(sim_out.ch_symb_log),'+')
  grid;

  figure(4)
  clf;
  plot(sim_out.ch_symb_log,'+')

  figure(5)
  clf;
  plot(sim_out.ratio_log)
  xlabel('time')
  ylabel('ratio for sync')
endfunction


function [freq_off_log EsNodBSet] = freq_off_est_test_curves
  EsNodBSet = [20 12 8];

  sim_in.frames    = 100;
  sim_in.foff      = -40;
  sim_in.dfoff     = 0;
  sim_in.fading_en = 1;
  freq_off_log = 1E6*ones(sim_in.frames, length(EsNodBSet) );
  sync_time_log = 1E6*ones(sim_in.frames, length(EsNodBSet) );

  for i=1:length(EsNodBSet)
    sim_in.EsNodB = EsNodBSet(i);
    printf("%f\n", sim_in.EsNodB);

    sim_out = freq_off_est_test(sim_in);
    freq_off_log(1:length(sim_out.freq_offset_log),i) = sim_out.freq_offset_log;
    sync_time_log(1:length(sim_out.sync_time_log),i) = sim_out.sync_time_log;
  end

  figure(1)
  clf
  hold on;
  for i=1:length(EsNodBSet)
    data = freq_off_log(find(freq_off_log(:,i) < 1E6),i);
    s = std(data);
    m = mean(data);
    stdbar = [m-s; m+s];
    plot(EsNodBSet(i), data, '+')
    plot([EsNodBSet(i) EsNodBSet(i)]+0.5, stdbar,'+-')
  end
  hold off

  axis([6 22 -25 25])
  if sim_in.fading_en
    title_str = sprintf('foff = %d Hz Fading', sim_in.foff);
  else
    title_str = sprintf('foff = %d Hz AWGN', sim_in.foff);
  end
  title(title_str);
  xlabel('Es/No (dB)')
  ylabel('freq offset error (Hz)');
 
  figure(2)
  clf
  hold on;
  for i=1:length(EsNodBSet)
    leg = sprintf("%d;%d dB;", i, EsNodBSet(i));
    plot(freq_off_log(find(freq_off_log(:,i) < 1E6),i),leg)
  end
  hold off;
  title(title_str);
  xlabel('test');
  ylabel('freq offset error (Hz)');
  legend('boxoff');

  figure(3)
  clf
  hold on;
  for i=1:length(EsNodBSet)
    data = sync_time_log(find(sync_time_log(:,i) < 1E6),i);
    if length(data) 
      s = std(data);
      m = mean(data);
      stdbar = [m-s; m+s];
      plot(EsNodBSet(i), data, '+')
      plot([EsNodBSet(i) EsNodBSet(i)]+0.5, stdbar,'+-')
    end
  end 
  hold off;
  axis([6 22 0 10])
  ylabel('sync time (frames)')
  xlabel('Es/No (dB)');
  title(title_str);

  figure(4)
  clf
  hold on;
  for i=1:length(EsNodBSet)
    leg = sprintf("%d;%d dB;", i, EsNodBSet(i));
    plot(sync_time_log(find(sync_time_log(:,i) < 1E6),i),leg)
  end
  hold off;
  title(title_str);
  xlabel('test');
  ylabel('sync time (frames)');
  legend('boxoff');

endfunction


% select on of these to run:

freq_off_est_test_single;
%freq_off_est_test_curves;

