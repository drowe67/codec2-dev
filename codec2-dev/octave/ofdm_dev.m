% ofdm_dev.m
% David Rowe April 2017
%
% Simulations used for development and testing of Rate Fs BPSK/QPSK
% OFDM modem.

ofdm_lib;
gp_interleaver;
ldpc;

#{
  TODO:
    [ ] run_sim neeeds to be refactored for coded operation at Nc=17 with UW
#}

function [sim_out rx states] = run_sim(sim_in)

  % set up core modem constants

  states = ofdm_init(sim_in.bps, sim_in.Rs, sim_in.Tcp, sim_in.Ns, sim_in.Nc);
  ofdm_load_const;
  Nbitspervocframe = 28;

  % simulation parameters and flags

  woffset = 2*pi*sim_in.foff_hz/Fs;
  dwoffset = 0;
  if isfield(sim_in, "dfoff_hz_per_sec")
    dwoffset = 2*pi*sim_in.dfoff_hz_per_sec/(Fs*Fs);
  end
  EbNodB  = sim_in.EbNodB;
  verbose = states.verbose = sim_in.verbose;
  hf_en   = sim_in.hf_en;

  timing_en = states.timing_en = sim_in.timing_en;
  states.foff_est_en = foff_est_en = sim_in.foff_est_en;
  states.phase_est_en = phase_est_en = sim_in.phase_est_en;
  if hf_en
    assert(phase_est_en == 1, "\nNo point running HF simulation without phase est!!\n");
  end
  if isfield(sim_in, "diversity_en")
    diversity_en = sim_in.diversity_en;
  else
    diversity_en = 0;
  end

  if verbose == 2
    printf("Rs:..........: %4.2f\n", Rs);
    printf("M:...........: %d\n", M);
    printf("Ncp:.........: %d\n", Ncp);
    printf("bps:.........: %d\n", bps);
    printf("Nbitsperframe: %d\n", Nbitsperframe);
    printf("Nrowsperframe: %d\n", Nrowsperframe);
    printf("Nsamperframe.: %d\n", Nsamperframe);
  end
  
  % Important to define run time in seconds so HF model will evolve the same way
  % for different pilot insertion rates.  So lets work backwards from approx
  % seconds in run to get Nbits, the total number of payload data bits

  Nrows = sim_in.Nsec*Rs;
  Nframes = floor((Nrows-1)/Ns);

  % if we are interleaving over multiple frames, adjust Nframes so we have an integer number
  % of interleaver frames in simulation

  interleave_en = 0;
  if isfield(sim_in, "interleave_frames") 
    interleave_frames = sim_in.interleave_frames;
    Nframes = interleave_frames*round(Nframes/interleave_frames);
    interleave_en = 1;
  end

  Nbits = Nframes * Nbitsperframe;    % number of payload data bits

  Nr = Nbits/(Nc*bps);                % Number of data rows to get Nbits total

  % double check if Nbits fit neatly into carriers

  assert(Nbits/(Nc*bps) == floor(Nbits/(Nc*bps)), "Nbits/(Nc*bps) must be an integer");

  Nrp = Nr + Nframes + 1;  % number of rows once pilots inserted
                           % extra row of pilots at end

  if verbose == 2
    printf("Nc...........: %d\n", Nc);
    printf("Ns...........: %d (step size for pilots, Ns-1 data symbols between pilots)\n", Ns);
    printf("Nr...........: %d\n", Nr);
    printf("Nbits........: %d\n", Nbits);
    printf("Nframes......: %d\n", Nframes);
    if interleave_en
      printf("Interleave fr: %d\n", interleave_frames);
    end
    printf("Nrp..........: %d (number of rows including pilots)\n", Nrp);
  end

  % Optional LPDC code -----------------------------------------------

  ldpc_en = states.ldpc_en = sim_in.ldpc_en;
  if sim_in.ldpc_en 
    assert(bps == 2, "Only QPSK supported for LDPC so far.....");
    HRA = sim_in.ldpc_code;
    [aNr aNc] = size(HRA);
    rate = states.rate = (aNc-aNr)/aNc;
    Ndatabitsperframe = Nbitsperframe;
    assert(aNc == Ndatabitsperframe, "Dude: Num cols of LDPC HRA must == Nbitsperframe");
    [H_rows, H_cols] = Mat2Hrows(HRA); 
    code_param.H_rows = H_rows; 
    code_param.H_cols = H_cols;
    code_param.P_matrix = [];
    code_param.data_bits_per_frame = length(code_param.H_cols) - length( code_param.P_matrix ); 
    code_param.code_bits_per_frame = aNc;
    assert(aNr == Ndatabitsperframe*rate);

    modulation = states.ldpc_modulation = 'QPSK';
    mapping = states.ldpc_mapping = 'gray';
    demod_type = states.ldpc_demod_type = 0;
    decoder_type = states.ldpc_decoder_type = 0;
    max_iterations = states.ldpc_max_iterations = 100;

    code_param.S_matrix = CreateConstellation( modulation, 4, mapping );

    states.code_param = code_param;
  elseif diversity_en
    rate = 0.5;    
  else
    rate = 1;
  end

  % set up HF model ---------------------------------------------------------------

  if hf_en

    % some typical values, or replace with user supplied

    dopplerSpreadHz = 1.0; path_delay_ms = 1;

    if isfield(sim_in, "dopplerSpreadHz") 
      dopplerSpreadHz = sim_in.dopplerSpreadHz;
    end
    if isfield(sim_in, "path_delay_ms") 
      path_delay_ms = sim_in.path_delay_ms;
    end
    path_delay_samples = path_delay_ms*Fs/1000;
    printf("Doppler Spread: %3.2f Hz Path Delay: %3.2f ms %d samples\n", dopplerSpreadHz, path_delay_ms, path_delay_samples);

    % generate same fading pattern for every run

    randn('seed',1);

    spread1 = doppler_spread(dopplerSpreadHz, Fs, Nrp*(M+Ncp)*1.1);
    spread2 = doppler_spread(dopplerSpreadHz, Fs, Nrp*(M+Ncp)*1.1);

    % sometimes doppler_spread() doesn't return exactly the number of samples we need
 
    assert(length(spread1) >= Nrp*(M+Ncp), "not enough doppler spreading samples");
    assert(length(spread2) >= Nrp*(M+Ncp), "not enough doppler spreading samples");
  end
 
  % ------------------------------------------------------------------
  % simulate for each Eb/No point 
  % ------------------------------------------------------------------

  for nn=1:length(EbNodB)
    rand('seed',1);
    randn('seed',1);

    EsNo = rate * bps * (10 .^ (EbNodB(nn)/10));
    variance = 1/(M*EsNo/2);

    Nsam = Nrp*(M+Ncp);

    % generate tx bits, optionaly LDPC encode, and modulate as QPSK symbols
    % note for reasons unknown LdpcEncode() returns garbage if we use > 0.5 rather than round()

    %tx_data_bits = round(rand(1,Nbits*rate));

    % std test frame so we can x-check
    
    tx_data_bits = create_ldpc_test_frame;
    
    tx_bits = []; tx_symbols = [];
    for f=1:Nframes
      st = (f-1)*Nbitsperframe*rate+1; en = st + Nbitsperframe*rate - 1;
      if ldpc_en
        codeword = LdpcEncode(tx_data_bits(st:en), code_param.H_rows, code_param.P_matrix);
      elseif diversity_en
        % Nc carriers, so Nc*bps bits/row coded, or Nc*bps*rate data bits that we repeat
        codeword = [];
        for rr=1:Ns-1
          st1 = st + (rr-1)*Nc*bps*rate; en1 = st1 + Nc*bps*rate - 1;
          codeword = [codeword tx_data_bits(st1:en1) tx_data_bits(st1:en1)];
        end
        assert(length(codeword) == Nbitsperframe);
      else
        % uncoded mode
        codeword = tx_data_bits;
        if isfield(sim_in, "uw_debug")
          codeword(states.uw_ind) = states.tx_uw;
        end
      end
      tx_bits = [tx_bits codeword];
      for b=1:2:Nbitsperframe
        tx_symbols = [tx_symbols qpsk_mod(codeword(b:b+1))];
      end
    end

    % optional interleaving over multiple frames

    if interleave_en
      for f=1:interleave_frames:Nframes
       st = (f-1)*Nbitsperframe/bps+1; en = st + Nbitsperframe*interleave_frames/bps - 1;
       tx_symbols(st:en) = gp_interleave(tx_symbols(st:en));
      end 
    end

    % OFDM transmitter

    tx = []; 
    for f=1:Nframes
      st = (f-1)*Nbitsperframe/bps+1; en = st + Nbitsperframe/bps - 1;
      tx = [tx ofdm_txframe(states, tx_symbols(st:en))];
    end
    
    % add extra row of pilots at end, to allow one frame simulations,
    % useful for development

    st = Nsamperframe*(Nframes-1)+1; en = st+Ncp+M-1;
    tx = [tx tx(st:en)];
    assert(length(tx) == Nsam);

    % channel simulation ---------------------------------------------------------------
    
    if isfield(sim_in, "sample_clock_offset_ppm") 
      % todo: this only works for large ppm like 500, runs out of memory
      %       for small ppm

      if sim_in.sample_clock_offset_ppm
        timebase = floor(abs(1E6/sim_in.sample_clock_offset_ppm));
        if sim_in.sample_clock_offset_ppm > 0
          tx = resample(tx, timebase+1, timebase);
        else
          tx = resample(tx, timebase, timebase+1);
        end

        % make sure length is correct for rest of simulation

        tx = [tx zeros(1,Nsam-length(tx))];
        tx = tx(1:Nsam);
      end
    end

    rx = tx;    

    if hf_en

      rx  = tx(1:Nsam) .* spread1(1:Nsam);
      rx += [zeros(1,path_delay_samples) tx(1:Nsam-path_delay_samples)] .* spread2(1:Nsam);

      % normalise rx power to same as tx

      nom_rx_pwr = 2/(Ns*(M*M)) + Nc/(M*M);
      rx_pwr = var(rx);
      rx *= sqrt(nom_rx_pwr/rx_pwr);
    end

    phase_offset = woffset*(1:Nsam) + 0.5*dwoffset*((1:Nsam).^2);
    rx = rx .* exp(j*phase_offset);
   
    if isfield(sim_in, "initial_noise_sams")
      rx = [zeros(1, sim_in.initial_noise_sams) rx];
      Nsam = length(rx);
    end
    
    noise = sqrt(variance)*(0.5*randn(1,Nsam) + j*0.5*randn(1,Nsam));
    snrdB = 10*log10(var(rx)/var(noise)) + 10*log10(8000) - 10*log10(3000);
    rx += noise;

    % interfering carrier

    % rx += 0.04*cos((1:length(rx))*states.w(10));
    
    % gain
    
    rx *= sim_in.gain;
    
    % some spare samples at end to avoid overflow as est windows may poke into the future a bit

    rx = [rx zeros(1,Nsamperframe)];
    
    % bunch of logs

    phase_est_pilot_log = [];
    delta_t_log = []; 
    timing_est_log = [];
    foff_est_hz_log = [];
    Nerrs_log = []; Nerrs_coded_log = [];
    rx_bits = []; rx_np = []; rx_amp = [];
    sig_var_log = []; noise_var_log = [];
    uw_errors_log = [];
    
    % reset some states for each EbNo simulation point

    states.sample_point = states.timing_est = 1;
    if timing_en == 0
      states.sample_point = Ncp;
    end
    states.nin = Nsamperframe;
    states.foff_est_hz = 0;    

    % for this simulation we "prime" buffer to allow one frame runs during development

    prx = 1;
    states.rxbuf(M+Ncp+2*Nsamperframe+1:Nrxbuf) = rx(prx:Nsamperframe+2*(M+Ncp));
    prx += Nsamperframe+2*(M+Ncp);

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

      [arx_bits states aphase_est_pilot_log arx_np arx_amp] = ofdm_demod(states, rxbuf_in);
     
      rx_bits = [rx_bits arx_bits]; rx_np = [rx_np arx_np]; rx_amp = [rx_amp arx_amp];

      % note: only supported in ldpc_en = 0 Nc=17 atm, see debug_false_sync()
      
      if isfield(sim_in, "uw_debug") 
        rx_uw = arx_bits(states.uw_ind);
        uw_errors = sum(xor(states.tx_uw,rx_uw));
        if verbose
          printf("f: %d uw_errors: %d\n", f, uw_errors);
        end
        uw_errors_log = [uw_errors_log uw_errors];
      end
      
      timing_est_log = [timing_est_log states.timing_est];
      delta_t_log = [delta_t_log states.delta_t];
      foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
      phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];
      sig_var_log = [sig_var_log states.sig_var];
      noise_var_log = [noise_var_log states.noise_var];
    end
    assert(length(rx_bits) == Nbits);

    % Optional de-interleave on rx QPSK symbols

    if interleave_en
      for f=1:interleave_frames:Nframes
       st = (f-1)*Nbitsperframe/bps+1; en = st + Nbitsperframe*interleave_frames/bps - 1;
       rx_np(st:en) = gp_deinterleave(rx_np(st:en));
       rx_amp(st:en) = gp_deinterleave(rx_amp(st:en));
      end 
    end

    % Calculate raw BER/PER stats, after pilots extracted.  As we may
    % have used interleaving, we qpsk_demod() here rather than using
    % rx_bits from ofdm_demod()

    rx_bits = zeros(1, Nbits);
    for s=1:Nbits/bps
      rx_bits(2*(s-1)+1:2*s) = qpsk_demod(rx_np(s));
    end

    errors = xor(tx_bits, rx_bits);
    Terrs = sum(errors);

    Tpackets = Tpacketerrs = 0;
    Nvocframes = floor(Nbits/Nbitspervocframe);
    for fv=1:Nvocframes
      st = (fv-1)*Nbitspervocframe + 1;
      en = st + Nbitspervocframe - 1;
      Nvocpacketerrs = sum(xor(tx_bits(st:en), rx_bits(st:en)));
      if Nvocpacketerrs
        Tpacketerrs++;
      end
      Tpackets++;
    end
 
    % Per-modem frame error log and optional LDPC/diversity error stats

    Terrs_coded = 0; Tpackets_coded = 0; Tpacketerrs_coded = 0; sim_out.error_positions = [];
    
    for f=1:Nframes
      st = (f-1)*Nbitsperframe+1; en = st + Nbitsperframe - 1;
      Nerrs_log(f) = sum(xor(tx_bits(st:en), rx_bits(st:en)));

      st = (f-1)*Nbitsperframe/bps + 1;
      en = st + Nbitsperframe/bps - 1;
      r = rx_np(st:en); fade = rx_amp(st:en);
      
      % optional LDPC decode
     
      if ldpc_en

        % scale based on amplitude ests

        mean_amp = states.mean_amp;
        rx_codeword = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, r/mean_amp, min(EsNo,30), fade/mean_amp);
      end

      % optional diversity demod

      if diversity_en
        rx_codeword = [];
        for rr=1:Ns-1
          for c=1:Nc/2
            s = (rr-1)*Nc + c;
            rx_codeword = [rx_codeword qpsk_demod(r(s)+r(s+Nc/2))];
          end
        end
        assert(length(rx_codeword) == Nbitsperframe*rate);
      end

      % running coded BER calcs

      if ldpc_en || diversity_en

        st = (f-1)*Nbitsperframe*rate + 1;
        en = st + Nbitsperframe*rate - 1;
        errors = xor(tx_data_bits(st:en), rx_codeword(1:Nbitsperframe*rate));
        Nerrs_coded = sum(errors);
        Nerrs_coded_log(f) = Nerrs_coded;
        Terrs_coded += Nerrs_coded;
        sim_out.error_positions = [sim_out.error_positions errors];

        % PER based on vocoder packet size, not sure it makes much
        % difference compared to using all bits in LDPC code for
        % packet

        atx_data_bits = tx_data_bits(st:en);
        Nvocframes = Nbitsperframe*rate/Nbitspervocframe;
        for fv=1:Nvocframes
          st = (fv-1)*Nbitspervocframe + 1;
          en = st + Nbitspervocframe - 1;
          Nvocpacketerrs = sum(xor(atx_data_bits(st:en), rx_codeword(st:en)));
          if Nvocpacketerrs
            Tpacketerrs_coded++;
          end
          Tpackets_coded++;
        end
      end
    end

    % print results of this simulation point to the console

    if verbose
      if ldpc_en || diversity_en
        printf("Coded EbNodB: % -4.1f BER: %5.4f Tbits: %5d Terrs: %5d PER: %5.4f Tpackets: %5d Tpacket_errs: %5d\n", 
                EbNodB(nn), Terrs_coded/(Nbits*rate), Nbits*rate, Terrs_coded, 
                Tpacketerrs_coded/Tpackets_coded, Tpackets_coded, Tpacketerrs_coded);
      end
      EbNodB_raw = EbNodB(nn) + 10*log10(rate);
      printf("Raw EbNodB..: % -4.1f BER: %5.4f Tbits: %5d Terrs: %5d PER: %5.4f Tpackets: %5d Tpacket_errs: %5d\n", 
             EbNodB_raw, Terrs/Nbits, Nbits, Terrs,
             Tpacketerrs/Tpackets, Tpackets, Tpacketerrs);
      EsNo = mean(sig_var_log)/mean(noise_var_log);
      %printf("Es/No est dB: % -4.1f\n", 10*log10(EsNo));
      sim_out.snrdB(nn) = snrdB;
      sim_out.snr_estdB(nn) = 10*log10(EsNo) + 10*log10(Nc*Rs/3000);
      
      % returns results for plotting curves

      if ldpc_en || diversity_en
        sim_out.ber(nn) = Terrs_coded/(Nbits*rate); 
        sim_out.per(nn) = Tpacketerrs_coded/Tpackets_coded; 
      else
        sim_out.ber(nn) = Terrs/Nbits; 
        sim_out.per(nn) = Tpacketerrs/Tpackets; 
      end
    end

    sim_out.uw_errors_log = uw_errors_log;
    
    % Optional plots, mostly used with run-single

    if verbose

      figure(1); clf; 
      plot(rx_np,'+');
      %axis([-2 2 -2 2]);
      title('Scatter');
      
      figure(2); clf;
      plot(phase_est_pilot_log,'g+', 'markersize', 5); 
      title('Phase est');
      axis([1 Nrp -pi pi]);  

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
      if ldpc_en
        subplot(211)
        stem(Nerrs_log/Nbitsperframe);
        title("Uncoded BER/frame");
        subplot(212)
        stem(Nerrs_coded_log/(Nbitsperframe*rate));
        title("Coded BER/frame");
      else
        title("BER/frame");
        stem(Nerrs_log/Nbitsperframe);
      end

      figure(6)
      Tx = abs(fft(rx(1:Nsam).*hanning(Nsam)'));
      Tx_dB = 20*log10(Tx);
      dF = Fs/Nsam;
      plot((1:Nsam)*dF, Tx_dB);
      mx = max(Tx_dB);
      axis([0 Fs/2 mx-60 mx])

      figure(7); clf;
      plot(10*log10(sig_var_log),'b;Es;');
      hold on;
      plot(10*log10(noise_var_log),'r;No;');
      snr_estdB = 10*log10(sig_var_log) - 10*log10(noise_var_log) + 10*log10(Nc*Rs/3000);
      snr_est_smoothed_dB = filter(0.1,[1 -0.9],snr_estdB);
      plot(snr_estdB,'g;SNR3k;');
      plot(snr_est_smoothed_dB,'c;SNR3k smooth;');
      
      hold off;
      title('Signal and Noise Power estimates');
      
#{
      if hf_en
        figure(4); clf; 
        subplot(211)
        plot(abs(spread1(1:Nsam)));
        %hold on; plot(abs(spread2(1:Nsam)),'g'); hold off;
        subplot(212)
        plot(angle(spread1(1:Nsam)));
        title('spread1 amp and phase');
      end
#}
      
#{
      % todo, work out a way to plot rate Fs hf model phase
      if sim_in.hf_en
        plot(angle(hf_model(:,2:Nc+1)));
      end
#}
    end

  end
endfunction


function run_single(EbNodB = 100, error_pattern_filename);
  Ts = 0.018; 
  sim_in.Tcp = 0.002; 
  sim_in.Rs = 1/Ts; sim_in.bps = 2; sim_in.Nc = 16; sim_in.Ns = 8;

  sim_in.Nsec = (sim_in.Ns+1)/sim_in.Rs;  % one frame, make sure sim_in.interleave_frames = 1
  sim_in.Nsec = 10;

  sim_in.EbNodB = 40;
  sim_in.verbose = 1;
  sim_in.hf_en = 0;
  sim_in.foff_hz = 0;
  sim_in.dfoff_hz_per_sec = 0.00;
  sim_in.sample_clock_offset_ppm = 0;
  sim_in.gain = 1;
  
  sim_in.timing_en = 0;
  sim_in.foff_est_en = 0;
  sim_in.phase_est_en = 1;

  load HRA_112_112.txt
  sim_in.ldpc_code = HRA_112_112;
  sim_in.ldpc_en = 1;

  sim_in.interleave_frames = 1;
  %sim_in.diversity_en = 1;

  sim_out = run_sim(sim_in);

  if nargin == 2
    fep = fopen(error_pattern_filename, "wb");
    fwrite(fep, sim_out.error_positions, "short");
    fclose(fep);
  end

end


% Plot BER and PER curves for AWGN and HF:
%
% i) BER/PER against Eb/No for various coding schemes
% ii) BER/PER against Eb/No showing Pilot/CP overhead
% iii) BER/PER against SNR, with pilot/CP overhead, comparing 700C

function run_curves

  % waveform

  Ts = 0.018; sim_in.Tcp = 0.002; 
  sim_in.Rs = 1/Ts; sim_in.bps = 2; sim_in.Nc = 16; sim_in.Ns = 8;

  pilot_overhead = (sim_in.Ns-1)/sim_in.Ns;
  cp_overhead = Ts/(Ts+sim_in.Tcp);
  overhead_dB = -10*log10(pilot_overhead*cp_overhead);

  % simulation parameters

  sim_in.verbose = 0;
  sim_in.foff_hz = 0;
  sim_in.gain = 1;
  
  sim_in.timing_en = 1;
  sim_in.foff_est_en = 1;
  sim_in.phase_est_en = 1;
  load HRA_112_112.txt
  sim_in.ldpc_code = HRA_112_112;
  sim_in.ldpc_en = 0;
  sim_in.hf_en = 0;

  sim_in.Nsec = 20;
  sim_in.EbNodB = 0:8;
  awgn_EbNodB = sim_in.EbNodB;

  awgn_theory = 0.5*erfc(sqrt(10.^(sim_in.EbNodB/10)));
  awgn = run_sim(sim_in);
  sim_in.ldpc_en = 1; awgn_ldpc = run_sim(sim_in);

  % Note for HF sim you really need >= 60 seconds (at Rs-50) to get sensible results c.f. theory

  sim_in.hf_en = 1; sim_in.ldpc_en = 0; 
  sim_in.Nsec = 60;
  sim_in.EbNodB = 4:2:14;

  EbNoLin = 10.^(sim_in.EbNodB/10);
  hf_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));

  hf = run_sim(sim_in);
  sim_in.diversity_en = 1; hf_diversity = run_sim(sim_in); sim_in.diversity_en = 0; 
  sim_in.ldpc_en = 1;  hf_ldpc = run_sim(sim_in);

  % try a few interleavers

  sim_in.interleave_frames = 1; hf_ldpc_1 = run_sim(sim_in);
  sim_in.interleave_frames = 8; hf_ldpc_8 = run_sim(sim_in);
  sim_in.interleave_frames = 16; hf_ldpc_16 = run_sim(sim_in);
  sim_in.interleave_frames = 32; hf_ldpc_32 = run_sim(sim_in);
  
  % Rate Fs modem BER curves of various coding schemes

  figure(1); clf;
  semilogy(awgn_EbNodB, awgn_theory,'b+-;AWGN theory;');
  hold on;
  semilogy(sim_in.EbNodB, hf_theory,'b+-;HF theory;');
  semilogy(awgn_EbNodB, awgn.ber,'r+-;AWGN;');
  semilogy(sim_in.EbNodB, hf.ber,'r+-;HF;');
  semilogy(sim_in.EbNodB, hf_diversity.ber,'ro-;HF diversity;');
  semilogy(awgn_EbNodB, awgn_ldpc.ber,'c+-;AWGN LDPC (224,112);');
  semilogy(sim_in.EbNodB, hf_ldpc.ber,'c+-;HF LDPC (224,112);');
  semilogy(sim_in.EbNodB, hf_ldpc_1.ber,'m+-;HF LDPC (224,112) interleave 1;');
  semilogy(sim_in.EbNodB, hf_ldpc_8.ber,'g+-;HF LDPC (224,112) interleave 8;');
  semilogy(sim_in.EbNodB, hf_ldpc_16.ber,'k+-;HF LDPC (224,112) interleave 16;');
  semilogy(sim_in.EbNodB, hf_ldpc_32.ber,'k+-;HF LDPC (224,112) interleave 32;');
  hold off;
  axis([0 14 1E-3 2E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
  legend("location", "southwest");
  title('Rate Fs modem BER for various FEC coding schemes');
  print('-deps', '-color', "ofdm_dev_ber_coding.eps")

  awgn_per_theory = 1 - (1-awgn_theory).^28;
  hf_per_theory = 1 - (1-hf_theory).^28;

  % Rate Fs modem PER curves of various coding schemes

  figure(2); clf;
  semilogy(awgn_EbNodB, awgn_per_theory,'b+-;AWGN theory;');
  hold on;
  semilogy(sim_in.EbNodB, hf_per_theory,'b+-;HF theory;');
  semilogy(awgn_EbNodB, awgn.per,'r+-;AWGN;');
  semilogy(sim_in.EbNodB, hf.per,'r+-;HF sim;');
  semilogy(sim_in.EbNodB, hf_diversity.per,'ro-;HF diversity;');
  semilogy(awgn_EbNodB, awgn_ldpc.per,'c+-;AWGN LDPC (224,112);');
  semilogy(sim_in.EbNodB, hf_ldpc.per,'c+-;HF LDPC (224,112);');
  semilogy(sim_in.EbNodB, hf_ldpc_1.per,'m+-;HF LDPC (224,112) interleave 1;');
  semilogy(sim_in.EbNodB, hf_ldpc_8.per,'g+-;HF LDPC (224,112) interleave 8;');
  semilogy(sim_in.EbNodB, hf_ldpc_16.per,'ko-;HF LDPC (224,112) interleave 16;');
  semilogy(sim_in.EbNodB, hf_ldpc_32.per,'k+-;HF LDPC (224,112) interleave 32;');
  hold off;
  axis([0 14 1E-2 1])
  xlabel('Eb/No (dB)');
  ylabel('PER');
  grid; grid minor on;
  legend('boxoff');
  legend("location", "southwest");
  title('Rate Fs modem PER for various FEC coding schemes');
  print('-deps', '-color', "ofdm_dev_per_coding.eps")

  % Rate Fs modem pilot/CP overhead BER curves

  figure(3); clf;
  semilogy(awgn_EbNodB, awgn_theory,'b+-;AWGN theory;');
  hold on;
  semilogy(sim_in.EbNodB, hf_theory,'b+-;HF theory;');
  semilogy(awgn_EbNodB+overhead_dB, awgn_theory,'g+-;AWGN lower bound pilot + CP;');
  semilogy(sim_in.EbNodB+overhead_dB, hf_theory,'g+-;HF lower bound pilot + CP;');
  semilogy(awgn_EbNodB+overhead_dB, awgn.ber,'r+-;AWGN sim;');
  semilogy(sim_in.EbNodB+overhead_dB, hf.ber,'r+-;HF sim;');
  hold off;
  axis([0 14 1E-3 2E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
  legend("location", "southwest");
  title('Rate Fs modem BER Pilot/Cyclic Prefix overhead');
  print('-deps', '-color', "ofdm_dev_ber_overhead.eps")

  % Rate Fs modem pilot/CP overhead PER curves

  figure(4); clf;
  semilogy(awgn_EbNodB, awgn_per_theory,'b+-;AWGN theory;','markersize', 10, 'linewidth', 2);
  hold on;
  semilogy(sim_in.EbNodB, hf_per_theory,'b+-;HF theory;','markersize', 10, 'linewidth', 2);
  semilogy(awgn_EbNodB+overhead_dB, awgn_per_theory,'g+-;AWGN lower bound pilot + CP;','markersize', 10, 'linewidth', 2);
  semilogy(sim_in.EbNodB+overhead_dB, hf_per_theory,'g+-;HF lower bound pilot + CP;','markersize', 10, 'linewidth', 2);
  semilogy(awgn_EbNodB+overhead_dB, awgn.per,'r+-;AWGN sim;','markersize', 10, 'linewidth', 2);
  semilogy(sim_in.EbNodB+overhead_dB, hf.per,'r+-;HF sim;','markersize', 10, 'linewidth', 2);
  hold off;
  axis([0 14 1E-2 1])
  xlabel('Eb/No (dB)');
  ylabel('PER');
  grid; grid minor on;
  legend('boxoff');
  legend("location", "southwest");
  title('Rate Fs modem PER Pilot/Cyclic Prefix overhead');
  print('-deps', '-color', "ofdm_dev_per_overhead.eps")

  % SNR including pilots, CP, and 700C est

  snr_awgn_theory = awgn_EbNodB + 10*log10(700/3000);
  snr_hf_theory   = sim_in.EbNodB + 10*log10(700/3000);
  snr_awgn        = snr_awgn_theory + overhead_dB;
  snr_hf          = sim_in.EbNodB + 10*log10(700/3000) + overhead_dB;

  % est 700C: 2/6 symbols are pilots, 1dB implementation loss

  snr_awgn_700c   = awgn_EbNodB + 10*log10(700/3000) + 10*log10(6/4) + 1;
  snr_hf_700c     = sim_in.EbNodB + 10*log10(700/3000) + 10*log10(6/4) + 1;

  figure(5); clf;
  semilogy(snr_awgn_theory, awgn_theory,'b+-;AWGN theory;','markersize', 10, 'linewidth', 2);
  hold on;
  semilogy(snr_awgn_700c, awgn_theory,'g+-;AWGN 700C;','markersize', 10, 'linewidth', 2);
  semilogy(snr_hf_700c, hf_diversity.ber,'go-;HF 700C;','markersize', 10, 'linewidth', 2);
  semilogy(snr_hf_theory, hf_theory,'b+-;HF theory;','markersize', 10, 'linewidth', 2);
  semilogy(snr_awgn, awgn_ldpc.ber,'c+-;AWGN LDPC (224,112);','markersize', 10, 'linewidth', 2);
  semilogy(snr_hf, hf_ldpc.ber,'c+-;HF LDPC (224,112);','markersize', 10, 'linewidth', 2);
  semilogy(snr_hf, hf_diversity.ber,'bo-;HF diversity;','markersize', 10, 'linewidth', 2);
  semilogy(snr_hf, hf_ldpc_16.ber,'k+-;HF LDPC (224,112) interleave 16;','markersize', 10, 'linewidth', 2);
  hold off;
  axis([-5 8 1E-3 2E-1])
  xlabel('SNR (3000Hz noise BW) (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
  legend("location", "southwest");
  title('Rate Fs modem BER versus SNR including pilot/CP overhead');
  print('-deps', '-color', "ofdm_dev_ber_snr.eps")

end


% Plot BER and PER curves for AWGN and HF with various estimators

function run_curves_estimators

  Nsec_awgn = 20;
  Nsec_hf = 60;

  % waveform

  Ts = 0.018; sim_in.Tcp = 0.002; 
  sim_in.Rs = 1/Ts; sim_in.bps = 2; sim_in.Nc = 16; sim_in.Ns = 8;

  % simulation parameters

  sim_in.verbose = 0; sim_in.foff_hz = 0; sim_in.ldpc_en = 0;

  sim_in.phase_est_en = 1;
  sim_in.timing_en = sim_in.foff_est_en = 0;

  % AWGN simulations

  sim_in.hf_en = 0; sim_in.Nsec = Nsec_awgn; sim_in.EbNodB = 0:2:6;
  sim_in.timing_en = sim_in.foff_est_en = 0;

  awgn_EbNodB = sim_in.EbNodB;
  awgn_theory = 0.5*erfc(sqrt(10.^(sim_in.EbNodB/10)));
  awgn = run_sim(sim_in);
  sim_in.timing_en = 1; awgn_timing = run_sim(sim_in);
  sim_in.foff_est_en = 1; awgn_foff_est = run_sim(sim_in);

  sim_in.dfoff_hz_per_sec = 0.02; awgn_dfoff = run_sim(sim_in);

  % HF simulations

  sim_in.hf_en = 1; sim_in.Nsec = Nsec_hf; sim_in.EbNodB = 4:2:8;
  sim_in.timing_en = sim_in.foff_est_en = 0;

  EbNoLin = 10.^(sim_in.EbNodB/10);
  hf_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));

  hf = run_sim(sim_in);
  sim_in.timing_en = 1; hf_timing = run_sim(sim_in);
  sim_in.timing_en = 0; sim_in.foff_est_en = 1; hf_foff_est = run_sim(sim_in);
  sim_in.timing_en = 1; hf_timing_foff_est = run_sim(sim_in);

  sim_in.dfoff_hz_per_sec = 0.02; hf_dfoff = run_sim(sim_in); sim_in.dfoff_hz_per_sec = 0.0;

  figure(1); clf;
  semilogy(awgn_EbNodB, awgn_theory,'b+-;AWGN theory;');
  hold on;
  semilogy(awgn_EbNodB, awgn.ber,'r+-;AWGN phase;');
  semilogy(awgn_EbNodB, awgn_timing.ber,'go-;AWGN phase+timing;');
  semilogy(awgn_EbNodB, awgn_foff_est.ber,'d+-;AWGN phase+timing+foff_est;');
  semilogy(awgn_EbNodB, awgn_dfoff.ber,'m+-;AWGN all + 0.02Hz/s drift;');
  semilogy(sim_in.EbNodB, hf_theory,'b+-;HF theory;');
  semilogy(sim_in.EbNodB, hf.ber,'r+-;HF phase;');
  semilogy(sim_in.EbNodB, hf_timing.ber,'go-;HF phase+timing;');
  semilogy(sim_in.EbNodB, hf_timing_foff_est.ber,'kd-;HF phase+foff_est;');
  semilogy(sim_in.EbNodB, hf_foff_est.ber,'c+-;HF phase+timing+foff_est;');
  semilogy(sim_in.EbNodB, hf_dfoff.ber,'m+-;HF + 0.02Hz/s drift;');
  hold off;
  axis([0 8 5E-3 1E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
  legend("location", "southeast");
  title('Rate Fs modem BER with estimators');
  print('-deps', '-color', "ofdm_dev_estimators.eps")

  % Simulate different HF path delays

  sim_in.hf_en = 1; sim_in.Nsec = Nsec_hf; sim_in.EbNodB = 4:2:8;
  sim_in.timing_en = sim_in.foff_est_en = 0;

  sim_in.path_delay_ms = 0; hf0 = run_sim(sim_in);
  sim_in.path_delay_ms = 0.5; hf500us = run_sim(sim_in);
  sim_in.path_delay_ms = 1; hf1ms = run_sim(sim_in);
  sim_in.path_delay_ms = 2; hf2ms = run_sim(sim_in);

  figure(2); clf;
  semilogy(sim_in.EbNodB, hf_theory,'b+-;HF theory;');
  hold on;
  semilogy(sim_in.EbNodB, hf0.ber,'r+-;HF phase 0 ms;');
  semilogy(sim_in.EbNodB, hf500us.ber,'g+-;HF phase 0.5 ms;');
  semilogy(sim_in.EbNodB, hf1ms.ber,'c+-;HF phase 1 ms;');
  semilogy(sim_in.EbNodB, hf2ms.ber,'k+-;HF phase 2 ms;');
  hold off;
  axis([3 9 1E-2 1E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
  legend("location", "southeast");
  title('Rate Fs modem BER across HF path delay');
  print('-deps', '-color', "ofdm_dev_estimators.eps")
end


% Plot SNR actual and estimated for various channels
%   Note no acquistion, as this can upset results, e.g. if
%   sync is lost.  We average SNR over entire run, in practice
%   there will be some sort of IIR averager.

function run_curves_snr

  Nsec_awgn = 30;
  Nsec_hf = 60;

  % waveform

  Ts = 0.018; sim_in.Tcp = 0.002; 
  sim_in.Rs = 1/Ts; sim_in.bps = 2; sim_in.Nc = 16; sim_in.Ns = 8;

  % simulation parameters

  sim_in.verbose = 1; sim_in.foff_hz = 0; sim_in.ldpc_en = 0; sim_in.gain = 1;
  sim_in.phase_est_en = sim_in.timing_en = sim_in.foff_est_en = 1;

  % AWGN simulation

  sim_in.EbNodB = 0:3:12;
  sim_in.hf_en = 0; sim_in.Nsec = Nsec_awgn;

  awgn = run_sim(sim_in);
  
  % HF simulations

  sim_in.hf_en = 1; sim_in.Nsec = Nsec_hf; sim_in.dopplerSpreadHz = 1; hf_fast = run_sim(sim_in);
  sim_in.dopplerSpreadHz = 0.5; hf_mid = run_sim(sim_in);
  sim_in.dopplerSpreadHz = 0.2; hf_slow = run_sim(sim_in);
  
  figure(1); clf;
  plot(sim_in.EbNodB, awgn.snrdB,'b+-;AWGN SNR;');
  hold on;
  plot(sim_in.EbNodB, awgn.snr_estdB   ,'g+-;AWGN SNR est;');
  plot(sim_in.EbNodB, hf_slow.snr_estdB,'k+-;HF 0.2 Hz SNR est;');
  plot(sim_in.EbNodB, hf_mid.snr_estdB ,'r+-;HF 0.5 Hz SNR est;');
  plot(sim_in.EbNodB, hf_fast.snr_estdB,'c+-;HF 1.0 Hz SNR est;');
  hold off;
  xlabel('Eb/No (dB)');
  ylabel('SNR (dB)');
  grid;
  legend('boxoff');
  title('SNR Actual and Estimated');
  legend("location", "northwest");
  print('-deps', '-color', "ofdm_dev_snr.eps")

end


% Run an acquisition test, returning vectors of estimation errors.
% Generates a vector of noise followed by continous rx signal to
% simulate signal starting.  Freq est will improve over time as we
% have an averaging statistic.  We then measure how long it takes to
% get good timing and freq estimates.

function [delta_ct delta_foff timing_mx_log] = acquisition_test(Ntests=10, EbNodB=100, foff_hz=0, hf_en=0, verbose=0)

  Ts = 0.018; 
  sim_in.Tcp = 0.002; 
  sim_in.Rs = 1/Ts; sim_in.bps = 2; sim_in.Nc = 17; sim_in.Ns = 8;

  sim_in.Nsec = (Ntests+1)*(sim_in.Ns+1)/sim_in.Rs;

  #{
    Notes:
      1) uncoded modem operating point, e.g -1dB for AWGN
      2) run_sim adds complex noise, when we take the real() below, this relects
         the -ve noise over to the +ve side, increasing the noise by 3dB.  In
         ofdm_tx.m and friends, we add real nosie that is correctly scaled so
         no problemo.  This means we need to increase the Eb/No below by 3dB
         to ensure the correct level of nosie ta the input of the timing_est.         
  #}
  
  sim_in.EbNodB = EbNodB + 3;
  sim_in.verbose = 0;
  sim_in.hf_en = hf_en;
  sim_in.foff_hz = foff_hz;
  sim_in.gain = 1;
  sim_in.timing_en = 1;
  sim_in.foff_est_en = 1;
  sim_in.phase_est_en = 1;
  sim_in.ldpc_en = 0;

  % optionally stick a bunch of noise in front of signal to confuse things
  
  sim_in.initial_noise_sams = 0;
  
  [sim_out rx states] = run_sim(sim_in);

  states.verbose = verbose;
  
  % set up acquistion 

  Nsamperframe = states.Nsamperframe; M = states.M; Ncp = states.Ncp;
  rate_fs_pilot_samples = states.rate_fs_pilot_samples;

  % test fine or acquisition over test signal

  delta_ct = []; delta_foff = []; timing_mx_log = []; foff_metric_log = [];

  
    % coarse acquiistion test.  We have no idea of timing or freq
    % offset for coarse we just use constant window shifts to simulate
    % a bunch of trials, this allows averaging of freq est
    % metric over time as we receive more and more frames

    st = 0.5*Nsamperframe; 
    en = 2.5*Nsamperframe - 1;    % note this gives Nsamperframe possibilities for coarse timing

    % actual known position of correct coarse timing

    ct_target = mod(sim_in.initial_noise_sams + Nsamperframe/2, Nsamperframe);

    i = 1;
    states.foff_metric = 0;
    for w=1:Nsamperframe:length(rx)-4*Nsamperframe
      [ct_est timing_valid timing_mx] = est_timing(states, real(rx(w+st:w+en)), rate_fs_pilot_samples);
      [foff_est states] = est_freq_offset(states, real(rx(w+st:w+en)), rate_fs_pilot_samples, ct_est);
      if states.verbose
        printf("i: %2d w: %5d ct_est: %4d foff_est: %5.1f timing_mx: %3.2f timing_vld: %d\n", i++, w, ct_est, foff_est, timing_mx, timing_valid);
      end

      % valid coarse timing ests are modulo Nsamperframe

      delta_ct = [delta_ct ct_est-ct_target];
      delta_foff = [delta_foff (foff_est-foff_hz)];
      timing_mx_log = [timing_mx_log; timing_mx];
      foff_metric_log = [foff_metric_log states.foff_metric];
    end

  if states.verbose > 1
    %printf("mean: %f std: %f\n", mean(delta_foff), std(delta_foff));
    figure(1); clf; plot(timing_mx_log,'+-');
    figure(2); clf; plot(delta_ct,'+-');
    figure(3); clf; plot(delta_foff,'+-');
    figure(4); clf; plot(foff_metric_log,'+');
    figure(5); clf; plot(real(rx))
  end
  
endfunction


#{

   Generates aquisistion statistics for AWGN and HF channels for
   continuous signals. Probability of acquistion is what matters,
   e.g. if it's 50% we can expect sync within 2 frames.

#}

function res = acquisition_histograms(fine_en = 0, foff, EbNoAWGN=-1, EbNoHF=3, verbose=1)
  Fs = 8000;
  Ntests = 100;

  % allowable tolerance for acquistion

  ftol_hz = 1.5;            % we can sync up on this
  ttol_samples = 0.002*Fs;  % 2ms, ie CP length

  % AWGN channel at uncoded Eb/No operating point

  [dct dfoff] = acquisition_test(Ntests, EbNoAWGN, foff, 0, fine_en);

  % Probability of acquistion is what matters, e.g. if it's 50% we can
  % expect sync within 2 frames

  PtAWGN = length(find (abs(dct) < ttol_samples))/length(dct);
  printf("AWGN P(time offset acq) = %3.2f\n", PtAWGN);
  if fine_en == 0
    PfAWGN = length(find (abs(dfoff) < ftol_hz))/length(dfoff);
    printf("AWGN P(freq offset acq) = %3.2f\n", PfAWGN);
  end

  if verbose
    figure(1); clf;
    hist(dct(find (abs(dct) < ttol_samples)))
    t = sprintf("Coarse Timing Error AWGN EbNo = %3.2f foff = %3.1f", EbNoAWGN, foff);
    title(t)
    if fine_en == 0
      figure(2)
      hist(dfoff(find(abs(dfoff) < 2*ftol_hz)))
      t = sprintf("Coarse Freq Error AWGN EbNo = %3.2f foff = %3.1f", EbNoAWGN, foff);
      title(t);
    end
 end

  % HF channel at uncoded operating point

  [dct dfoff] = acquisition_test(Ntests, EbNoHF, foff, 1, fine_en);

  PtHF = length(find (abs(dct) < ttol_samples))/length(dct);
  printf("HF P(time offset acq) = %3.2f\n", PtHF);
  if fine_en == 0
    PfHF = length(find (abs(dfoff) < ftol_hz))/length(dfoff)
    printf("HF P(freq offset acq) = %3.2f\n", PfHF);
  end

  if verbose
    figure(3); clf;
    hist(dct(find (abs(dct) < ttol_samples)))
    t = sprintf("Coarse Timing Error HF EbNo = %3.2f foff = %3.1f", EbNoHF, foff);
    title(t)
    if fine_en == 0
      figure(4)
      hist(dfoff(find(abs(dfoff) < 2*ftol_hz)))
      t = sprintf("Coarse Freq Error HF EbNo = %3.2f foff = %3.1f", EbNoHF, foff);
      title(t);
    end
  end
  
  res = [PtAWGN PfAWGN PtHF PfHF];
endfunction


% plot some curves of Acquisition probability against EbNo and freq offset

function acquistion_curves

  EbNo = [-1 2 5 8];
  %foff = [-20 -15 -10 -5 0 5 10 15 20];
  foff = [-15 -5 0 5 15];
  cc = ['b' 'g' 'k' 'c' 'm'];
  
  figure(1); clf; hold on; title('P(timing) AWGN'); xlabel('Eb/No dB'); legend('location', 'southeast');
  figure(2); clf; hold on; title('P(freq) AWGN'); xlabel('Eb/No dB'); legend('location', 'southeast');
  figure(3); clf; hold on; title('P(timing) HF'); xlabel('Eb/No dB'); legend('location', 'southeast');
  figure(4); clf; hold on; title('P(freq) HF'); xlabel('Eb/No dB'); legend('location', 'southeast');
  
  for f = 1:length(foff)
    afoff = foff(f);
    res_log = [];
    for e = 1:length(EbNo)
      aEbNo = EbNo(e);
      res = zeros(1,4);
      res = acquisition_histograms(fine_en = 0, afoff, aEbNo, aEbNo+4, verbose = 0);
      res_log = [res_log; res];
    end
    figure(1); l = sprintf('%c+-;%3.1f Hz;', cc(f), afoff); plot(EbNo, res_log(:,1), l);
    figure(2); l = sprintf('%c+-;%3.1f Hz;', cc(f), afoff); plot(EbNo, res_log(:,3), l);
    figure(3); l = sprintf('%c+-;%3.1f Hz;', cc(f), afoff); plot(EbNo+4, res_log(:,2), l);
    figure(4); l = sprintf('%c+-;%3.1f Hz;', cc(f), afoff); plot(EbNo+4, res_log(:,4), l);
  end

  figure(1); print('-deps', '-color', "ofdm_dev_acq_curves_time_awgn.eps")
  figure(2); print('-deps', '-color', "ofdm_dev_acq_curves_freq_awgn.eps")
  figure(3); print('-deps', '-color', "ofdm_dev_acq_curves_time_hf.eps")
  figure(4); print('-deps', '-color', "ofdm_dev_acq_curves_freq_hf.eps")
endfunction


% Used to develop sync state machine - in particular a metric to show
% we are out of sync, or have sync with a bad freq offset est, or have
% lost modem signal

function sync_metrics(x_axis = 'EbNo')
  Fs      = 8000;
  Ntests  = 4;
  f_offHz = [-25:25];
  EbNodB  = [-10 0 3 6 10 20];
  %f_offHz = [-5:5:5];
  %EbNodB  = [-10 0 10];
  cc = ['b' 'g' 'k' 'c' 'm' 'b'];
  pt = ['+' '+' '+' '+' '+' 'o'];
  
  mean_mx1_log = mean_dfoff_log = [];
  for f = 1:length(f_offHz)
    af_offHz = f_offHz(f);
    mean_mx1_row = mean_dfoff_row = [];
    for e = 1:length(EbNodB)
      aEbNodB = EbNodB(e);
      [dct dfoff timing_mx_log] = acquisition_test(Ntests, aEbNodB, af_offHz);
      mean_mx1 = mean(timing_mx_log(:,1));
      printf("f_offHz: %5.2f EbNodB: % 6.2f mx1: %3.2f\n", af_offHz, aEbNodB, mean_mx1);
      mean_mx1_row = [mean_mx1_row mean_mx1];
      mean_dfoff_row = [mean_dfoff_row mean(dfoff)];
    end
    mean_mx1_log = [mean_mx1_log; mean_mx1_row];
    mean_dfoff_log = [mean_dfoff_log; mean_dfoff_row];
  end

  figure(1); clf; hold on; grid;
  if strcmp(x_axis,'EbNo')
    for f = 1:length(f_offHz)
      if f == 2, hold on, end;
      leg1 = sprintf("b+-;mx1 %4.1f Hz;", f_offHz(f));
      plot(EbNodB, mean_mx1_log(f,:), leg1)
    end
    hold off;
    xlabel('Eb/No (dB)');
    ylabel('Coefficient')
    title('Pilot Correlation Metric against Eb/No for different Freq Offsets');
    legend("location", "northwest"); legend("boxoff");
    axis([min(EbNodB) max(EbNodB) 0 1.2])
    print('-deps', '-color', "ofdm_dev_pilot_correlation_ebno.eps")
  end

  if strcmp(x_axis,'freq')
    % x axis is freq

    for e = length(EbNodB):-1:1
      leg1 = sprintf("%c%c-;mx1 %3.0f dB;", cc(e), pt(e), EbNodB(e));
      plot(f_offHz, mean_mx1_log(:,e), leg1)
    end
    hold off;
    xlabel('freq offset (Hz)');
    ylabel('Coefficient')
    title('Pilot Correlation Metric against Freq Offset for different Eb/No dB');
    legend("location", "northwest"); legend("boxoff");
    axis([min(f_offHz) max(f_offHz) 0 1])
    print('-deps', '-color', "ofdm_dev_pilot_correlation_freq.eps")

    mean_dfoff_log
 
    figure(2); clf;
    for e = 1:length(EbNodB)
      if e == 2, hold on, end;
      leg1 = sprintf("+-;mx1 %3.0f dB;", EbNodB(e));
      plot(f_offHz, mean_dfoff_log(:,e), leg1)
    end
    hold off;
    xlabel('freq offset (Hz)');
    ylabel('Mean Freq Est Error')
    title('Freq Est Error against Freq Offset for different Eb/No dB');
    axis([min(f_offHz) max(f_offHz) -5 5])
  end
  
endfunction


% during development it was discovered demod could obtain a flase sync with no UW
% errors at +/- 7Hz, approx the frame rate.  This function is used to explore that

function debug_false_sync(EbNodB = 100)
  Ts = 0.018; 
  sim_in.Tcp = 0.002; 
  sim_in.Rs = 1/Ts; sim_in.bps = 2; sim_in.Nc = 17; sim_in.Ns = 8;

  sim_in.Nsec = (sim_in.Ns+1)/sim_in.Rs;  % one frame, make sure sim_in.interleave_frames = 1
  sim_in.Nsec = 1;

  sim_in.EbNodB = 40;
  sim_in.verbose = 0;
  sim_in.hf_en = 0;
  sim_in.foff_hz = 0;
  sim_in.dfoff_hz_per_sec = 0.00;
  sim_in.sample_clock_offset_ppm = 0;
  sim_in.gain = 1;
  
  sim_in.timing_en = 0;
  sim_in.foff_est_en = 0;
  sim_in.phase_est_en = 1;

  load HRA_112_112.txt
  sim_in.ldpc_code = HRA_112_112;
  sim_in.ldpc_en = 0;

  %sim_in.interleave_frames = 1;
  %sim_in.diversity_en = 1;

  sim_in.uw_debug = 1;

  foff = -25:0.5:25;
  for f=1:length(foff)
    sim_in.foff_hz = foff(f);
    sim_out = run_sim(sim_in);
    min_uw_errors(f) = min(sim_out.uw_errors_log(2:end-1));
    printf("f: %4.2f %2d\n", foff(f), min_uw_errors(f));
  end
  figure(1); clf;
  plot(foff, min_uw_errors);
  title('UW errors versus freq offset');
  xlabel('Freq Offset (Hz)');
end


% Using an input raw file, plot frame by frame metric information,
% used to debug false syncs

function metric_fbf(fn, Nsec)
  Ts = 0.018; 
  states = ofdm_init(bps=2, Rs=1/Ts, Tcp=0.002, Ns=8, Nc=17);
  ofdm_load_const;
  states.verbose = 2;

  % factor of 2 as input is a real valued signal

  Ascale = states.amp_scale/2; 
  f = fopen(fn,"rb"); rx = fread(f,Inf,"short")'/Ascale; fclose(f);
  if (nargin == 2) && (length(rx) > Nsec*Fs)
    rx = rx(1:Nsec*Fs);
  end
  Nsam = length(rx);
  %bpf_coeff = make_ofdm_bpf(write_c_header_file=0);
  %rx = filter(bpf_coeff,1,rx);
  
  st = 0.5*Nsamperframe; 
  en = 2.5*Nsamperframe - 1;    % note this gives Nsamperframe possibilities for coarse timing

  i = 1; w_log = timing_mx_log = av_level_log = [];
  states.foff_metric = 0;
  for w=1:Nsamperframe:length(rx)-4*Nsamperframe
    printf("%3d %5d", i,w);
    if i == 30
      states.verbose = 3;
    else
      states.verbose = 2;
    end  
    [ct_est timing_valid timing_mx av_level] = est_timing(states, real(rx(w+st:w+en)), states.rate_fs_pilot_samples);
    i++;
    w_log = [w_log w];
    timing_mx_log = [timing_mx_log timing_mx];
    av_level_log = [av_level_log av_level];
  end

  figure(2); clf;
  mx = max(abs(rx)); 
  subplot(211); plot(rx); axis([0 Nsam -mx mx]);
  subplot(212); hold on;
  plot(w_log,timing_mx_log,'b+-;timing mx;');
  plot(w_log,av_level_log,'g+-;av level;');
  hold off;
endfunction


% ---------------------------------------------------------
% choose simulation to run here 
% ---------------------------------------------------------

format;
more off;

init_cml('~/cml/');

%run_single(10);
%run_curves
%run_curves_estimators
%acquisition_histograms(fin_en=0, foff_hz=-15, EbNoAWGN=-1, EbNoHF=3)
%acquisition_test(Ntests=10, EbNodB=-1, foff_hz=0, hf_en=0, verbose=1);
%sync_metrics('freq')
%run_curves_snr
%acquisition_dev(Ntests=10, EbNodB=100, foff_hz=0)
%acquistion_curves
%debug_false_sync
metric_fbf("~/Desktop/false_sync.wav")
