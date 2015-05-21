% test_cohpsk.m
% David Rowe Oct 2014
%

% Bunch of simulations to test coherent FDM QPSK with pilot based
% coherent detection, DSSS (diversity), and rate 1/2 LDPC

% TODO
%   [X] Nc carriers, 588 bit frames
%   [X] FEC
%   [X] pilot insertion and removal
%   [ ] delay on parity carriers
%   [X] pilot based phase est
%   [ ] uncoded and coded frame sync
%   [X] timing estimation, RN filtering, carrier FDM
%   [ ] this file getting too big - refactor
%   [ ] robust coarse timing
%   [ ] Np knob on GUI
%   [ ] experimental tuning on EnNo_, fading[], to optimise LDPC dec perf
 
% reqd to make sure we can repeat tests exactly

rand('state',1); 
randn('state',1);
graphics_toolkit ("gnuplot");

cohpsk;

function test_curves
     
  sim_in = standard_init();
  sim_in.do_write_pilot_file = 0;

  % single test point ---------------------------------------

  sim_in.verbose          = 1;
  sim_in.plot_scatter     = 1;

  sim_in.Esvec            = 20; 
  sim_in.framesize        = 32;
  sim_in.Ntrials          = 400;
  sim_in.Rs               = 50;
  sim_in.Nc               = 4;
  sim_in.Np               = 2;
  sim_in.Ns               = 4;
  sim_in.Nd               = 1;
  sim_in.modulation       = 'qpsk';
  sim_in.ldpc_code_rate   = 1;
  sim_in.ldpc_code        = 0;
  sim_in.coh_en           = 1;

  sim_in.hf_sim           = 1;
  sim_in.hf_mag_only      = 0;
  sim_in.f_off            = 0;
  sim_in.div_timeshift    = 1;

  sim_qpsk                = ber_test(sim_in);

  % AWGN curves ----------------------------------------------------

  sim_in.Ntrials          = 400;
  sim_in.hf_sim           = 0;
  sim_in.plot_scatter     = 0;
  sim_in.Esvec            = 5:10; 

  Ebvec = sim_in.Esvec - 10*log10(2);
  BER_theory = 0.5*erfc(sqrt(10.^(Ebvec/10)));
  
  sim_in.modulation       = 'dqpsk';
  sim_dqpsk               = ber_test(sim_in, 'dqpsk');

  sim_in.modulation       = 'qpsk';
  sim_qpsk_pilot          = ber_test(sim_in, 'qpsk');

  % HF curves ----------------------------------------------------

  sim_in.Ntrials          = 400;
  sim_in.hf_sim           = 1;
  sim_in.plot_scatter     = 0;
  sim_in.Esvec            = 5:20; 
  sim_in.modulation       = 'dqpsk';
  sim_dqpsk_hf            = ber_test(sim_in, 'dqpsk');

  sim_in.modulation       = 'qpsk';
  sim_in.coh_en           = 0;
  sim_in.hf_mag_only      = 1;
  sim_qpsk_hf             = ber_test(sim_in, 'qpsk');

  sim_in.coh_en           = 1;
  sim_in.hf_mag_only      = 0;
  sim_qpsk_pilot_hf       = ber_test(sim_in, 'qpsk');

  sim_in.Nd               = 2;
  sim_in.div_timeshift    = 1;
  sim_qpsk_pilot_hf_div   = ber_test(sim_in, 'qpsk');

  sim_in.div_timeshift    = sim_in.Rs;
  sim_qpsk_pilot_hf_div2  = ber_test(sim_in, 'qpsk');

  sim_in.modulation       = 'qpsk';
  sim_in.coh_en           = 0;
  sim_in.hf_mag_only      = 1;
  sim_in.div_timeshift    = 1;
  sim_qpsk_hf_div         = ber_test(sim_in, 'qpsk');

  % plot results ---------------------------------------------------

  figure(1); 
  clf;
  semilogy(Ebvec, BER_theory,'r;QPSK theory;')
  hold on;

  semilogy(sim_dqpsk.Ebvec, sim_dqpsk.BERvec,'c;DQPSK AWGN;')
  semilogy(sim_qpsk_pilot.Ebvec, sim_qpsk_pilot.BERvec,'b;QPSK pilot AWGN;')

  semilogy(sim_qpsk_hf.Ebvec, sim_qpsk_hf.BERvec,'r;QPSK HF ideal;')
  semilogy(sim_qpsk_hf_div.Ebvec, sim_qpsk_hf_div.BERvec,'r;QPSK HF ideal div;')
  semilogy(sim_qpsk_pilot_hf.Ebvec, sim_qpsk_pilot_hf.BERvec,'b;QPSK pilot HF;')
  semilogy(sim_qpsk_pilot_hf_div.Ebvec, sim_qpsk_pilot_hf_div.BERvec,'g;QPSK pilot Nd=2 HF;')
  semilogy(sim_qpsk_pilot_hf_div2.Ebvec, sim_qpsk_pilot_hf_div2.BERvec,'g;QPSK pilot Nd=2 Rs HF;')
  semilogy(sim_dqpsk_hf.Ebvec, sim_dqpsk_hf.BERvec,'c;DQPSK HF;')

  hold off;

  xlabel('Eb/N0')
  ylabel('BER')
  grid("minor")
  axis([min(sim_qpsk_hf.Ebvec) max(sim_qpsk_hf.Ebvec) 1E-3 1])
  legend("boxoff");
endfunction


function test_single

  sim_in = standard_init();

  sim_in.verbose          = 1;
  sim_in.plot_scatter     = 1;

  sim_in.framesize        = 32;
  sim_in.Nc               = 4;
  sim_in.Rs               = 50;
  sim_in.Ns               = 4;
  sim_in.Np               = 2;
  sim_in.Nd               = 2;
  sim_in.ldpc_code_rate   = 1;
  sim_in.ldpc_code        = 0;

  sim_in.Ntrials          = 400;
  sim_in.Esvec            = 12; 
  sim_in.hf_sim           = 1;
  sim_in.hf_mag_only      = 0;
  sim_in.modulation       = 'qpsk';
  sim_in.coh_en           = 1;
  sim_in.f_off            = 0;

  sim_in.div_timeshift    = 1;

  %sim_in.modulation      = 'dqpsk';

  sim_qpsk_hf             = ber_test(sim_in);

  %fep=fopen("errors_450.bin","wb"); fwrite(fep, sim_qpsk_hf.ldpc_errors_log, "short"); fclose(fep);
endfunction


% Rate Fs test funcs -----------------------------------------------------------

function rate_Fs_tx(tx_filename)
  sim_in = standard_init();

  sim_in.verbose          = 1;
  sim_in.plot_scatter     = 1;

  sim_in.framesize        = 160;
  sim_in.Nc               = 4;
  sim_in.Rs               = 50;
  sim_in.Ns               = 4;
  sim_in.Np               = 2;
  sim_in.Nchip            = 1;
  sim_in.ldpc_code_rate   = 1;
  sim_in.ldpc_code        = 0;

  sim_in.Ntrials          = 20;
  sim_in.Esvec            = 7; 
  sim_in.hf_sim           = 1;
  sim_in.hf_mag_only      = 0;
  sim_in.modulation       = 'qpsk';

  sim_in = symbol_rate_init(sim_in);

  prev_sym_tx             = sim_in.prev_sym_tx;
  code_param              = sim_in.code_param;
  tx_bits_buf             = sim_in.tx_bits_buf;
  framesize               = sim_in.framesize;
  rate                    = sim_in.ldpc_code_rate;
  Ntrials                 = sim_in.Ntrials;
  Rs                      = sim_in.Rs;
  Fs                      = sim_in.Fs;
  Nc                      = sim_in.Nc;

  M = Fs/Rs;

  EsNodB = sim_in.Esvec(1);
  EsNo = 10^(EsNodB/10);
 
  rx_symb_log = []; av_tx_pwr = [];

  rn_coeff = gen_rn_coeffs(0.5, 1/Fs, Rs, 6, M);
  tx_symb_buf = [];

  tx_bits = round(rand(1,framesize*rate));                       

  for nn=1:Ntrials+2

    [tx_symb tx_bits_out prev_sym_tx] = symbol_rate_tx(sim_in, tx_bits, code_param, prev_sym_tx);
    tx_bits_buf(1:framesize) = tx_bits_buf(framesize+1:2*framesize);
    tx_bits_buf(framesize+1:2*framesize) = tx_bits_out;
    tx_symb_buf = [tx_symb_buf; tx_symb];
  end
 
  % zero pad and tx filter

  [m n] = size(tx_symb_buf);
  zp = [];
  for i=1:m
    zrow = M*tx_symb_buf(i,:);
    zp = [zp; zrow; zeros(M-1,Nc)];
  end

  for c=1:Nc
    tx_filt(:,c) = filter(rn_coeff, 1, zp(:,c));
  end

  % upconvert to real IF and save to disk

  [m n] = size(tx_filt);
  tx_fdm = zeros(1,m);
  Fc = 1500;
  for c=1:Nc
    freq(c) = exp(j*2*pi*(Fc - c*Rs*1.5)/Fs);
  end
  phase_tx = ones(1,Nc);
  %phase_tx = exp(j*2*pi*(0:Nc)/(Nc+1));

  for c=1:Nc
    for i=1:m
      phase_tx(c) = phase_tx(c) * freq(c);
      tx_fdm(i) = tx_fdm(i) + tx_filt(i,c)*phase_tx(c);
    end
  end

  tx_fdm = real(tx_fdm);

  %tx_fdm = compress(tx_fdm, 0.4);
  %tx_fdm = sign(tx_fdm) .* (abs(tx_fdm) .^ 0.4); 
  %hpa_clip = max(abs(tx_fdm))*0.8
  %tx_fdm(find(abs(tx_fdm) > hpa_clip)) = hpa_clip;

  papr = max(tx_fdm.*conj(tx_fdm)) / mean(tx_fdm.*conj(tx_fdm));
  papr_dB = 10*log10(papr);
  printf("PAPR: %4.2f dB\n", papr_dB);

  Ascale = 2000;
  figure(1);
  clf;
  plot(Ascale*tx_fdm(1:8000))

  ftx=fopen(tx_filename,"wb"); fwrite(ftx, Ascale*real(tx_fdm), "short"); fclose(ftx);

endfunction


function rate_Fs_rx(rx_filename)
  sim_in = standard_init();

  sim_in.verbose          = 1;
  sim_in.plot_scatter     = 1;

  sim_in.framesize        = 160;
  sim_in.Nc               = 4;
  sim_in.Rs               = 50;
  sim_in.Ns               = 4;
  sim_in.Np               = 4;
  sim_in.Nchip            = 1;
  sim_in.ldpc_code_rate   = 1;
  sim_in.ldpc_code        = 0;

  sim_in.Ntrials          = 10;
  sim_in.Esvec            = 40; 
  sim_in.hf_sim           = 1;
  sim_in.hf_mag_only      = 0;
  sim_in.modulation       = 'qpsk';

  sim_in = symbol_rate_init(sim_in);

  prev_sym_tx             = sim_in.prev_sym_tx;
  prev_sym_rx             = sim_in.prev_sym_rx;
  code_param              = sim_in.code_param;
  tx_bits_buf             = sim_in.tx_bits_buf;
  framesize               = sim_in.framesize;
  rate                    = sim_in.ldpc_code_rate;
  Ntrials                 = sim_in.Ntrials;
  Rs                      = sim_in.Rs;
  Fs                      = sim_in.Fs;
  Nc                      = sim_in.Nc;
  Nsymbrowpilot           = sim_in.Nsymbrowpilot;
  pilot                   = sim_in.pilot;
  Ns                      = sim_in.Ns;
  Npilotsframe            = sim_in.Npilotsframe;

  M = Fs/Rs;

  EsNodB = sim_in.Esvec(1);
  EsNo = 10^(EsNodB/10);
 
  phi_log = []; amp_log = []; EsNo__log = [];
  rx_symb_log = []; av_tx_pwr = [];
  Terrs = Tbits = 0;
  errors_log = []; Nerrs_log = []; 
  ldpc_Nerrs_log = []; ldpc_errors_log = [];
  Ferrsldpc = Terrsldpc = Tbitsldpc = 0;

  rn_coeff = gen_rn_coeffs(0.5, 1/Fs, Rs, 6, M);

  tx_bits = round(rand(1,framesize*rate));                       

  % read from disk

  Ascale = 2000;
  frx=fopen(rx_filename,"rb"); rx_fdm = fread(frx, "short")/Ascale; fclose(frx);

  rx_fdm=sqrt(2)*rx_fdm;

  if 0 
    % optionally add AWGN noise for testing calibration of Es//No measurement

    snr_3000Hz_dB = -9;
    snr_4000Hz_lin = 10 ^ (snr_3000Hz_dB/10);
    snr_4000Hz_lin *= (3000/4000);
    variance = var(rx_fdm)/snr_4000Hz_lin;
    rx_fdm += sqrt(variance)*randn(length(rx_fdm),1);
  end

  figure(1)
  plot(rx_fdm(800:1200));

  % freq offset estimation

  printf("Freq offset and coarse timing est...\n");
  [f_max max_s_Fs] = test_freq_off_est(rx_filename, 1,5*6400);
  f_max = 0; max_s_Fs = 4;
  max_s = floor(max_s_Fs/M + 6);

  printf("Downconverting...\n");

  [m n] = size(rx_fdm);
  rx_symb = zeros(m,Nc);
  Fc = 1500;
  for c=1:Nc
    freq(c) = exp(-j*2*pi*(f_max + Fc - c*Rs*1.5)/Fs);
  end
  phase_rx = ones(1,Nc);
  rx_bb = zeros(m,Nc);

  for c=1:Nc
    for i=1:m
      phase_rx(c) = phase_rx(c) * freq(c);
      rx_bb(i,c) = rx_fdm(i)*phase_rx(c);
    end
  end

  printf("Filtering...\n");
  for c=1:Nc
    rx_filt(:,c) = filter(rn_coeff, 1, rx_bb(:,c));
  end

  %subplot(211);
  %plot(real(rx_filt(1:10*M,9)));
  %subplot(212);
  %plot(imag(rx_filt(1:10*M,9)));

  % Fine timing estimation and decimation to symbol rate Rs. Break rx
  % signal into ft sample blocks.  If clock offset is 1000ppm,
  % that's one more/less sample over Ft samples at Fs=8000 Hz.

  printf("Fine timing estimation....\n");
  ft = M*10;
  [nsam m] = size(rx_filt);
  rx_symb_buf = []; rx_timing_log = [];
  
  for st=1:ft:floor(nsam/ft - 1)*ft
    % fine timing and decimation

    env = zeros(ft,1);
    for c=1:Nc
      env = env + abs(rx_filt(st:st+ft-1,c));
    end

    % The envelope has a frequency component at the symbol rate.  The
    % phase of this frequency component indicates the timing.  So work out
    % single DFT at frequency 2*pi/M

    x = exp(-j*2*pi*(0:ft-1)/M) * env;
  
    norm_rx_timing = angle(x)/(2*pi);
    %norm_rx_timing = -0.4;
    if norm_rx_timing < 0
      rx_timing = -floor(norm_rx_timing*M+0.5) + M;
    else
      rx_timing = -floor(norm_rx_timing*M+0.5) + 2*M;
    end
    rx_timing_log = [rx_timing_log norm_rx_timing];

    % printf("%d %d\n", st+rx_timing, st+rx_timing+ft-1);
    rx_symb_buf = [rx_symb_buf; rx_filt(st+rx_timing:M:st+rx_timing+ft-1,:)];
  end
  
  figure(2)
  clf;
  plot(rx_timing_log)
  axis([1 length(rx_timing_log) -0.5 0.5 ])
  title('fine timing')
    
  printf("Symbol rate demodulation....\n");
  phase_off = 1;
  Ntrials = floor((nsam/M)/Nsymbrowpilot) - 2;
  %max_s = 6;

  for nn=1:Ntrials

    s_ch = rx_symb_buf((nn-1)*Nsymbrowpilot+max_s:nn*Nsymbrowpilot+max_s-1,:);
    [rx_symb rx_bits rx_symb_linear amp_linear amp_ phi_ EsNo_ prev_sym_rx sim_in] = symbol_rate_rx(sim_in, s_ch, prev_sym_rx);
        
    rx_symb_log = [rx_symb_log rx_symb_linear];
    phi_log = [phi_log; phi_];
    amp_log = [amp_log; amp_];

    if nn > 1
      EsNo__log = [EsNo__log EsNo_];

      % Measure BER

      error_positions = xor(rx_bits(1:framesize*rate), tx_bits(1:framesize*rate));
      Nerrs = sum(error_positions);
      Terrs += Nerrs;
      Tbits += framesize*rate;
      errors_log = [errors_log error_positions];
      Nerrs_log = [Nerrs_log Nerrs];

      if sim_in.ldpc_code
        % LDPC decode
            
        detected_data = ldpc_dec(code_param, sim_in.max_iterations, sim_in.demod_type, sim_in.decoder_type, ...
                                 rx_symb_linear, min(100,EsNo_), amp_linear);
        error_positions = xor(detected_data(1:framesize*rate), tx_bits(1:framesize*rate) );
        Nerrs = sum(error_positions);
        ldpc_Nerrs_log = [ldpc_Nerrs_log Nerrs];
        ldpc_errors_log = [ldpc_errors_log error_positions];
        if Nerrs
          Ferrsldpc++;
        end
        Terrsldpc += Nerrs;
        Tbitsldpc += framesize*rate;
      end
    end
  end

  EsNo_av = mean(10*log10(EsNo__log));
  printf("EsNo est (dB): %3.1f SNR est: %3.2f Terrs: %d BER %4.2f QPSK BER theory %4.2f av_tx_pwr: %3.2f", 
         EsNo_av, mean(EsNo_to_SNR(10*log10(EsNo__log))),
         Terrs,
         Terrs/Tbits, 0.5*erfc(sqrt(EsNo/2)), av_tx_pwr);
  if sim_in.ldpc_code
    printf("\n LDPC: Terrs: %d BER: %4.2f Ferrs: %d FER: %4.2f\n", 
           Terrsldpc, Terrsldpc/Tbitsldpc, Ferrsldpc, Ferrsldpc/Ntrials);
  end

  figure(3);
  clf;
  scat = rx_symb_log .* exp(j*pi/4);
  scat = scat((framesize):length(scat));
  plot(real(scat), imag(scat),'+');
  title('Scatter plot');
        
  figure(4)
  clf
  subplot(211)
  stem(Nerrs_log)
  axis([0 Ntrials+1 0 max(Nerrs_log)+1])
  title('Uncoded Errors/Frame');
  if sim_in.ldpc_code
    subplot(212)
    stem(ldpc_Nerrs_log)
    axis([0 Ntrials+1 0 max(ldpc_Nerrs_log)+1])
    title('Coded Errors/Frame');
  end

  figure(5)
  clf
  [m1 n1] = size(phi_log);
  phi_x = [];
  phi_x_counter = 1;
  p = Ns;
  for r=1:m1
    if p == Ns
      phi_x_counter++;
      p = 0;
    end
    p++;
    phi_x = [phi_x phi_x_counter++];        
  end
          
  subplot(211)
  plot(phi_x, phi_log(:,1),'r+;Estimated HF channel phase;')
  ylabel('Phase (rads)');
  legend('boxoff');

  subplot(212)
  plot(phi_x, amp_log(:,1),'r+;Estimated HF channel amp;')
  ylabel('Amplitude');
  xlabel('Time (symbols)');
  legend('boxoff');

  figure(6);
  clf
  plot(EsNo_to_SNR(10*log10(EsNo__log)));
  title('SNR est (dB)')

  fep=fopen("errors_450.bin","wb"); fwrite(fep, ldpc_errors_log, "short"); fclose(fep);

endfunction




function snr = EsNo_to_SNR(EsNo)
  snr = interp1([20 11.8 9.0 6.7 4.9 3.3 -10], [ 3 3 0 -3 -6 -9 -9], EsNo);
endfunction


function gen_test_bits()
  sim_in = standard_init();
  framesize = 32*10;
  tx_bits = round(rand(1,framesize));
  test_bits_coh_file(tx_bits);
endfunction




% Start simulations ---------------------------------------

more off;
%close all;
%test_curves();
test_single();
%rate_Fs_tx("tx_zero.raw");
%rate_Fs_tx("tx.raw");
%rate_Fs_rx("tx_-4dB.wav")
%rate_Fs_rx("tx.raw")
%gen_test_bits();

