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
    [ ] compute SNR and PAPR
    [ ] SSB bandpass filtering
    [ ] way to simulate aquisition and demod
    [ ] testframe based, maybe repeat every 10 seconds
        + work out which pattern we match to sync up
    [ ] acquisition curves
        + plot error versus freq and timing offset
        + plot pro acquisition versus freq offset, timing and freq sep and together
        + plot total acquist prob at various SNRs ... maybe mean number of frames to sync?
    [ ] replace genie EsNo est used for ldpc dec
    [ ] clean up text 
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
    assert(aNc == Nbitsperframe, "Dude: Num cols of LDPC HRA must == Nbitsperframe");
    [H_rows, H_cols] = Mat2Hrows(HRA); 
    code_param.H_rows = H_rows; 
    code_param.H_cols = H_cols;
    code_param.P_matrix = [];
    code_param.data_bits_per_frame = length(code_param.H_cols) - length( code_param.P_matrix ); 
    code_param.code_bits_per_frame = aNc;
    assert(aNr == Nbitsperframe*rate);

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

    tx_data_bits = round(rand(1,Nbits*rate));

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
        codeword = tx_data_bits(st:en);
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
   
    noise = sqrt(variance)*(0.5*randn(1,Nsam) + j*0.5*randn(1,Nsam));
    snrdB = 10*log10(var(rx)/var(noise)) + 10*log10(8000) - 10*log10(3000);
    rx += noise;

    % some spare samples at end to avoid overflow as est windows may poke into the future a bit

    rx = [rx zeros(1,Nsamperframe)];
    
    % bunch of logs

    phase_est_pilot_log = [];
    delta_t_log = []; 
    timing_est_log = [];
    foff_est_hz_log = [];
    Nerrs_log = []; Nerrs_coded_log = [];
    rx_bits = []; rx_np = []; rx_amp = [];

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
      timing_est_log = [timing_est_log states.timing_est];
      delta_t_log = [delta_t_log states.delta_t];
      foff_est_hz_log = [foff_est_hz_log states.foff_est_hz];
      phase_est_pilot_log = [phase_est_pilot_log; aphase_est_pilot_log];
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
        % note we put ceiling on EsNo as decoder misbehaives with high EsNo

        rx_codeword = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, r, min(EsNo,30), fade);
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

    if ldpc_en || diversity_en
      printf("Coded EbNodB: % -4.1f BER: %5.4f Tbits: %5d Terrs: %5d PER: %5.4f Tpackets: %5d Tpacket_errs: %5d\n", 
              EbNodB(nn), Terrs_coded/(Nbits*rate), Nbits*rate, Terrs_coded, 
              Tpacketerrs_coded/Tpackets_coded, Tpackets_coded, Tpacketerrs_coded);
    end
    EbNodB_raw = EbNodB(nn) + 10*log10(rate);
    printf("Raw EbNodB..: % -4.1f BER: %5.4f Tbits: %5d Terrs: %5d PER: %5.4f Tpackets: %5d Tpacket_errs: %5d\n", 
           EbNodB_raw, Terrs/Nbits, Nbits, Terrs,
           Tpacketerrs/Tpackets, Tpackets, Tpacketerrs);

    % returns results for plotting curves

    if ldpc_en || diversity_en
      sim_out.ber(nn) = Terrs_coded/(Nbits*rate); 
      sim_out.per(nn) = Tpacketerrs_coded/Tpackets_coded; 
    else
      sim_out.ber(nn) = Terrs/Nbits; 
      sim_out.per(nn) = Tpacketerrs/Tpackets; 
    end

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
      Tx = abs(fft(tx(1:Nsam).*hanning(Nsam)'));
      Tx_dB = 20*log10(Tx);
      dF = Fs/Nsam;
      plot((1:Nsam)*dF, Tx_dB);
      mx = max(Tx_dB);
      axis([0 Fs/2 mx-60 mx])

     
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
  %sim_in.Nsec = 1;

  sim_in.EbNodB = EbNodB;
  sim_in.verbose = 2;
  sim_in.hf_en = 0;
  sim_in.foff_hz = 0;
  sim_in.dfoff_hz_per_sec = 0.00;
  sim_in.sample_clock_offset_ppm = 0;

  sim_in.timing_en = 0;
  sim_in.foff_est_en = 0;
  sim_in.phase_est_en = 0;

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


% Run an acquisition test, returning vectors of estimation errors

function [delta_ct delta_foff] = acquisition_test(Ntests=10, EbNodB=100, foff_hz=0, hf_en=0, fine_en=0)

  Ts = 0.018; 
  sim_in.Tcp = 0.002; 
  sim_in.Rs = 1/Ts; sim_in.bps = 2; sim_in.Nc = 16; sim_in.Ns = 8;

  sim_in.Nsec = (Ntests+1)*(sim_in.Ns+1)/sim_in.Rs;

  sim_in.EbNodB = EbNodB;
  sim_in.verbose = 0;
  sim_in.hf_en = hf_en;
  sim_in.foff_hz = foff_hz; 
  sim_in.timing_en = 1;
  sim_in.foff_est_en = 1;
  sim_in.phase_est_en = 1;
  sim_in.ldpc_en = 0;
  
  [sim_out rx states] = run_sim(sim_in);
  
  states.verbose = 0;
  
  % set up acquistion 

  Nsamperframe = states.Nsamperframe; M = states.M; Ncp = states.Ncp;
  rate_fs_pilot_samples = states.rate_fs_pilot_samples;

  % test fine or acquisition over test signal
  #{
    fine: - start with coarse timing instant
          - on each frame est timing a few samples about that point
          - update timing instant

    corr: - where is best plcase to sample
          - just before end of symbol?
          - how long should sequence be?
          - add extra?
          - aim for last possible moment?
          - man I hope IL isn't too big.....
  #}

  delta_ct = []; delta_foff = [];

  % a fine simulation is a bit like what ofsd_demod() does, just searches a few samples
  % either side of current coarse est
  
  if fine_en

    window_width = 5;                   % search +/-2 samples from current timing instant
    timing_instant = Nsamperframe+1;    % start at correct instant for AWGN
                                        % start at second frame so we can search -2 ... +2

    while timing_instant < (length(rx) - (Nsamperframe + length(rate_fs_pilot_samples) + window_width))
      st = timing_instant - ceil(window_width/2); 
      en = st + Nsamperframe-1 + length(rate_fs_pilot_samples) + window_width-1;
      [ft_est foff_est] = coarse_sync(states, rx(st:en), rate_fs_pilot_samples);
      printf("ft_est: %d timing_instant %d %d\n", ft_est, timing_instant, mod(timing_instant, Nsamperframe));
      timing_instant += Nsamperframe + ft_est - ceil(window_width/2);
      delta_t = [delta_ft ft_est - ceil(window_width/2)];
    end
  else
    % for coarse simulation we just use constant window shifts

    st = 0.5*Nsamperframe; 
    en = 2.5*Nsamperframe - 1;    % note this gives Nsamperframe possibilities for coarse timing
    ct_target = Nsamperframe/2;   % actual known position of correct coarse timing

    for w=1:Nsamperframe:length(rx)-4*Nsamperframe
    %for w=1:M+Ncp:length(rx)-4*Nsamperframe
      [ct_est foff_est] = coarse_sync(states, rx(w+st:w+en), rate_fs_pilot_samples);
      if states.verbose
        printf("w: %d ct_est: %4d foff_est: %3.1f\n", w, ct_est, foff_est);
      end

      % valid coarse timing ests are modulo Nsamperframe

      delta_ct = [delta_ct ct_est-ct_target];
      delta_foff = [delta_foff (foff_est-foff_hz)];
    end
  end

endfunction


#{
   Run some tests for various acquisition conditions. Probability of
   acquistion is what matters, e.g. if it's 50% we can expect sync
   within 2 frames

  Results on 17 Mar 2018:
  
  foff Hz  Channel Eb/No  P(t)  P(f)
  -20      AWGN    -1     0.99  0.42
  -20      HF       3     0.94  0.43
  +20      AWGN    -1     1.00  0.37
  +20      HF       3     0.93  0.24
  -5       AWGN    -1     1.00  0.56
  -5       HF       3     0.98  0.50
  +5       AWGN    -1     1.00  0.54
  +5       HF       3     0.98  0.39
   0       AWGN    10     1.00  0.98
   0       HF      10     1.00  0.65

   -> Suggests we will sync up in 2-3 frames which is pretty cool.  Would be good
      to have freq est about as reliable as timing est.....
#}

function acquisition_histograms(fine_en = 0, foff)
  Fs = 8000;
  Ntests = 100;

  % allowable tolerance for acquistion

  ftol_hz = 2.0;            % fine freq can track this out
  ttol_samples = 0.002*Fs;  % 2ms, ie CP length

  % AWGN channel operating point

  [dct dfoff] = acquisition_test(Ntests, -1, foff, 0, fine_en);

  % Probability of acquistion is what matters, e.g. if it's 50% we can
  % expect sync within 2 frames

  printf("AWGN P(time offset acq) = %3.2f\n", length(find (abs(dct) < ttol_samples))/length(dct));
  if fine_en == 0
    printf("AWGN P(freq offset acq) = %3.2f\n", length(find (abs(dfoff) < ftol_hz))/length(dfoff));
  end

  figure(1)
  hist(dct(find (abs(dct) < ttol_samples)))
  title('Coarse Timing error AWGN')
  if fine_en == 0
    figure(2)
    hist(dfoff)
  title('Coarse Freq error AWGN')
  end

  % HF channel operating point

  [dct dfoff] = acquisition_test(Ntests, -3, foff, 1, fine_en);

  printf("HF P(time offset acq) = %3.2f\n", length(find (abs(dct) < ttol_samples))/length(dct));
  if fine_en == 0
    printf("HF P(freq offset acq) = %3.2f\n", length(find (abs(dfoff) < ftol_hz))/length(dfoff));
  end

  figure(3)
  hist(dct(find (abs(dct) < ttol_samples)))
  title('Coarse Timing error HF')
  if fine_en == 0
    figure(4)
    hist(dfoff)
    title('Coarse Freq error HF')
  end

endfunction


% ---------------------------------------------------------
% choose simulation to run here 
% ---------------------------------------------------------

format;
more off;

init_cml('/home/david/Desktop/cml/');

run_single 
%run_curves
%run_curves_estimators
%acquisition_histograms(0, 0)
%acquisition_test(10, 4, 5)
