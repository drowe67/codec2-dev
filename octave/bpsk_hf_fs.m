% bpsk_hf_fs.m
% David Rowe Mar 2017
%
% Rate Fs BPSK simulation, development of bpsk_hf_rs.m

#{
 TODO:
   [X] strip back experimental stuff to just features we need
   [X] ZOH/integrator
   [X] OFDM up and down conversion
   [X] rate Fs HF model and HF results
   [X] add QPSK
   [X] add CP
   [ ] adjust waveform parameters for real world
   [ ] Nsec run time take into account CP
   [ ] handle border carriers
       [ ] start with phantom carriers 
           + but unhappy with 1800Hz bandwidth
       [ ] also try interpolation or just single row
   [ ] compute SNR and PAPR
   [ ] acquisition & freq offset estimation
   [ ] SSB bandpass filtering
#}

1;

% Gray coded QPSK modulation function

function symbol = qpsk_mod(two_bits)
    two_bits_decimal = sum(two_bits .* [2 1]); 
    switch(two_bits_decimal)
        case (0) symbol =  1;
        case (1) symbol =  j;
        case (2) symbol = -j;
        case (3) symbol = -1;
    endswitch
endfunction


% Gray coded QPSK demodulation function

function two_bits = qpsk_demod(symbol)
    bit0 = real(symbol*exp(j*pi/4)) < 0;
    bit1 = imag(symbol*exp(j*pi/4)) < 0;
    two_bits = [bit1 bit0];
endfunction


function sim_out = run_sim(sim_in)
  Rs = 62.5;
  Fs = 8000;
  M  = Fs/Rs;
  Tcp = 0.002; Ncp = Tcp*Fs;
  foffset = 0;
  woffset = 2*pi*foffset/Fs;

  bps = sim_in.bps;
  EbNodB  = sim_in.EbNodB;
  verbose = sim_in.verbose;
  hf_en   = sim_in.hf_en;

  Ns = sim_in.Ns;          % step size for pilots
  Nc = sim_in.Nc;          % Number of cols, aka number of carriers

  Nbitsperframe = (Ns-1)*Nc*bps;
  Nrowsperframe = Nbitsperframe/(Nc*bps);
  if verbose
    printf("Rs:.........: %4.2f\n", Rs);
    printf("M:..........: %d\n", M);
    printf("bps:.........: %d\n", bps);
    printf("Nbitsperframe: %d\n", Nbitsperframe);
    printf("Nrowsperframe: %d\n", Nrowsperframe);
  end

  % Important to define run time in seconds so HF model will evolve the same way
  % for different pilot insertion rates.  So lets work backwards from approx
  % seconds in run to get Nbits, the total number of payload data bits

  % frame has Ns-1 data symbols between pilots, e.g. for Ns=3: 
  %
  % PPP
  % DDD
  % DDD
  % PPP

  Nrows = sim_in.Nsec*Rs
  Nframes = floor((Nrows-1)/Ns);
  Nbits = Nframes * Nbitsperframe;    % number of payload data bits

  Nr = Nbits/(Nc*bps);                % Number of data rows to get Nbits total

  if verbose
    printf("Nc.....: %d\n", Nc);
    printf("Ns.....: %d (step size for pilots, Ns-1 data symbols between pilots)\n", Ns);
    printf("Nr.....: %d\n", Nr);
    printf("Nbits..: %d\n", Nbits);
  end

  % double check if Nbits fit neatly into carriers

  assert(Nbits/(Nc*bps) == floor(Nbits/(Nc*bps)), "Nbits/(Nc*bps) must be an integer");
 
  printf("Nframes: %d\n", Nframes);

  Nrp = Nr + Nframes + 1;  % number of rows once pilots inserted
                           % extra row of pilots at end
  printf("Nrp....: %d (number of rows including pilots)\n", Nrp);

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

    randn('seed',1);
    spread1 = doppler_spread(dopplerSpreadHz, Fs, (sim_in.Nsec*(M+Ncp)/M+0.2)*Fs);
    spread2 = doppler_spread(dopplerSpreadHz, Fs, (sim_in.Nsec*(M+Ncp)/M+0.2)*Fs);

    % sometimes doppler_spread() doesn't return exactly the number of samples we need
 
    assert(length(spread1) >= Nrp*M, "not enough doppler spreading samples");
    assert(length(spread2) >= Nrp*M, "not enough doppler spreading samples");

    hf_gain = 1.0/sqrt(var(spread1)+var(spread2));
    % printf("nsymb: %d lspread1: %d\n", nsymb, length(spread1));
  end
 
  % simulate for each Eb/No point ------------------------------------

  for nn=1:length(EbNodB)
    rand('seed',1);
    randn('seed',1);

    EbNo = bps * (10 .^ (EbNodB(nn)/10));
    variance = 1/(M*EbNo/2);

    % generate tx bits

    tx_bits = rand(1,Nbits) > 0.5;

    % map to symbols in linear array

    if bps == 1
      tx_sym_lin = 2*tx_bits - 1;
    end
    if bps == 2
      for s=1:Nbits/bps
        tx_sym_lin(s) = qpsk_mod(tx_bits(2*(s-1)+1:2*s));
      end
    end

    % place symbols in multi-carrier frame with pilots and boundary carriers

    tx_sym = []; s = 1;
    for f=1:Nframes
      aframe = zeros(Nrowsperframe,Nc+2);
      aframe(1,:) = 1;
      for r=1:Nrowsperframe
        arowofsymbols = tx_sym_lin(s:s+Nc-1);
        s += Nc;
        aframe(r+1,2:Nc+1) = arowofsymbols;
      end
      tx_sym = [tx_sym; aframe];
    end      
    tx_sym = [tx_sym; ones(1,Nc+2)];  % final row of pilots
    [nr nc] = size(tx_sym);
    assert(nr == Nrp);

    % OFDM up conversion and upsampling to rate Fs ---------------------

    w = (0:Nc+1)*2*pi*Rs/Fs;
    W = zeros(Nc+2,M);
    for c=1:Nc+2
      W(c,:) = exp(j*w(c)*(0:M-1));
    end

    Nsam = Nrp*(M+Ncp);
    tx = [];

    % OFDM upconvert symbol by symbol so we can add CP

    for r=1:Nrp
      asymbol = tx_sym(r,:) * W/M;
      asymbol_cp = [asymbol(M-Ncp+1:M) asymbol];
      tx = [tx asymbol_cp];
    end
    assert(length(tx) == Nsam);

    % channel simulation ---------------------------------------------------------------

    rx = tx;
    
    if hf_en
      %rx  =  [zeros(1,path_delay_samples) tx(1:Nsam-path_delay_samples)];
      rx  = hf_gain * tx(1:Nsam) .* spread1(1:Nsam);
      rx  += hf_gain * [zeros(1,path_delay_samples) tx(1:Nsam-path_delay_samples)] .* spread2(1:Nsam);
    end

    rx = rx .* exp(j*woffset*(1:Nsam));

    noise = sqrt(variance)*(0.5*randn(1,Nsam) + j*0.5*randn(1,Nsam));
    rx += noise;
    
    % downconvert, downsample and integrate using OFDM.  Start integrating just after CP 
    % when ISI has died down

    rx_sym = zeros(Nrp, Nc+2);
    for r=1:Nrp
      st = (r-1)*(M+Ncp)+Ncp+1; en = st + M - 1;
      %printf("st: %d en: %d\n", st, en);
      for c=1:Nc+2
        acarrier = rx(st:en) .* conj(W(c,:));
        rx_sym(r,c) = sum(acarrier);
      end
    end
    
    % pilot based phase est, we use known tx symbols as pilots ----------
 
    phase_est_pilot_log = 10*ones(Nrp,Nc+2);
    phase_est_stripped_log = 10*ones(Nrp,Nc+2);
    phase_est_log = 10*ones(Nrp,Nc+2);
    for c=2:Nc+1
      for r=1:Ns:Nrp-Ns

        % estimate phase using average of 6 pilots in a rect 2D window centred
        % on this carrier
        % PPP
        % DDD
        % DDD
        % PPP
          
        cr = c-1:c+1;
        aphase_est_pilot_rect = sum(rx_sym(r,cr)*tx_sym(r,cr)') + sum(rx_sym(r+Ns,cr)*tx_sym(r+Ns,cr)');

        % use next step of pilots in past and future

        if r > Ns+1
          aphase_est_pilot_rect += sum(rx_sym(r-Ns,cr)*tx_sym(r-Ns,cr)');
        end
        if r < Nrp - 2*Ns
          aphase_est_pilot_rect += sum(rx_sym(r+2*Ns,cr)*tx_sym(r+2*Ns,cr)');
        end

        % correct phase offset using phase estimate

        for rr=r+1:r+Ns-1
          aphase_est_pilot = angle(aphase_est_pilot_rect);
          phase_est_pilot_log(rr,c) = aphase_est_pilot;
          rx_corr(rr,c) = rx_sym(rr,c) * exp(-j*aphase_est_pilot);
        end

      end % r=1:Ns:Nrp-Ns
    end % c=2:Nc+1

    % remove pilots to give us just data symbols and demodulate

    rx_bits = []; rx_np = [];
    for r=1:Nrp
      if mod(r-1,Ns) != 0
        arowofsymbols = rx_corr(r,2:Nc+1);
        rx_np = [rx_np arowofsymbols];
        if bps == 1
          arowofbits = real(arowofsymbols) > 0;
        end
        if bps == 2
          arowofbits = zeros(1,Nc);
          for c=1:Nc
            arowofbits((c-1)*2+1:c*2) = qpsk_demod(arowofsymbols(c));
          end
        end
        rx_bits = [rx_bits arowofbits];
      end
    end
    assert(length(rx_bits) == Nbits);


    % calculate BER stats as a block, after pilots extracted

    errors = xor(tx_bits, rx_bits);
    Nerrs = sum(errors);

    printf("EbNodB: %3.2f BER: %4.3f Nbits: %d Nerrs: %d\n", EbNodB(nn), Nerrs/Nbits, Nbits, Nerrs);

    if verbose
      figure(1)
      plot(real(tx))
      figure(2)
      Tx = abs(fft(tx.*hanning(Nsam)'));
      Tx_dB = 20*log10(Tx);
      dF = Fs/Nsam;
      plot((1:Nsam)*dF, Tx_dB);
      mx = max(Tx_dB);
      axis([0 Fs/2 mx-60 mx])
     
      figure(3); clf; 
      plot(rx_np,'+');
      axis([-2 2 -2 2]);

      
      if hf_en
        figure(4); clf; 
        subplot(211)
        plot(abs(spread1(1:Nsam)));
        %hold on; plot(abs(spread2(1:Nsam)),'g'); hold off;
        subplot(212)
        plot(angle(spread1(1:Nsam)));
      end
      

      figure(5); clf;
      plot(phase_est_log(:,2:Nc+1),'+', 'markersize', 10); 
      hold on; 
      plot(phase_est_pilot_log(:,2:Nc+1),'g+', 'markersize', 5); 

#{
      if sim_in.hf_en
        plot(angle(hf_model(:,2:Nc+1)));
      end
#}

      axis([1 Nrp -pi pi]);
    end

    sim_out.ber(nn) = sum(Nerrs)/Nbits; 
    sim_out.pilot_overhead = 10*log10(Ns/(Ns-1));
  end
endfunction


% Plot BER against Eb/No curves for AWGN and HF

% Target operating point Eb/No for HF is 6dB, as this is where our rate 1/2
% LDPC code gives good results (10% PER, 1% BER).  However this means
% the Eb/No at the input is 10*log(1/2) or 3dB less, so we need to
% make sure phase est works at Eb/No = 6 - 3 = 3dB
%
% For AWGN target is 2dB so -1dB op point.

function run_curves
  sim_in.bps = 2; sim_in.Nc = 8; sim_in.Ns = 7; sim_in.verbose = 0;

  sim_in.hf_en = 0;
  sim_in.Nsec = 30;
  sim_in.EbNodB = -3:5;
  awgn_EbNodB = sim_in.EbNodB;

  awgn_theory = 0.5*erfc(sqrt(10.^(sim_in.EbNodB/10)));
  awgn = run_sim(sim_in);

  sim_in.hf_en = 1;
  sim_in.Nsec = 120;
  sim_in.EbNodB = 1:8;

  EbNoLin = 10.^(sim_in.EbNodB/10);
  hf_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));

  hf = run_sim(sim_in);

  figure(4); clf;
  semilogy(awgn_EbNodB, awgn_theory,'b+-;AWGN theory;');
  hold on;
  semilogy(sim_in.EbNodB, hf_theory,'g+-;HF theory;');
  semilogy(awgn_EbNodB, awgn.ber,'r+-;AWGN sim;');
  semilogy(sim_in.EbNodB, hf.ber,'c+-;HF sim;');
  hold off;
  axis([-3 8 1E-2 2E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
end


function run_single
  sim_in.bps = 2;
  sim_in.Nc = 8;
  sim_in.Ns = 7;
  % sim_in.Nsec = (sim_in.Ns+1)/62.5; % one frame
  sim_in.Nsec = 120;
  sim_in.EbNodB = 3;
  sim_in.verbose = 1;
  sim_in.hf_en = 1;
  sim_in.path_delay_ms = 1;

  run_sim(sim_in);
end


format;
more off;

%run_single
run_curves




