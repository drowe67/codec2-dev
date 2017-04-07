% bpsk_hf_fs.m
% David Rowe Mar 2017
%
% Rate Fs BPSK simulation, development of bpsk_hf_rs.m

#{
 TODO:
   [X] strip back experimental stuff to just features we need
   [X] ZOH/integrator
   [X] OFDM up and down conversion
   [ ] rate Fs HF model and HF results
   [ ] handle boarder carriers
       [ ] start with phantom carriers 
           + but unhappy with 1800Hz bandwidth
       [ ] also try interpolation or just single row
   [ ] compute SNR and PAPR
   [ ] acquisition & freq offset estimation
   [ ] SSB bandpass filtering
#}

1;

function sim_out = run_sim(sim_in)
  Rs = 100;
  Fs = 8000;
  M  = Fs/Rs;

  EbNodB  = sim_in.EbNodB;
  verbose = sim_in.verbose;
  hf_en   = sim_in.hf_en;

  Ns = sim_in.Ns;          % step size for pilots
  Nc = sim_in.Nc;          % Number of cols, aka number of carriers

  Nbitsperframe = (Ns-1)*Nc;
  printf("Nbitsperframe: %d\n", Nbitsperframe);
  Nrowsperframe = Nbitsperframe/Nc;
  printf("Nrowsperframe: %d\n", Nrowsperframe);

  % Important to define run time in seconds so HF model will evolve the same way
  % for different pilot insertion rates.  So lets work backwards from approx
  % seconds in run to get Nbits, the total number of payload data bits

  % frame has Ns-1 data symbols between pilots, e.g. for Ns=3: 
  %
  % PPP
  % DDD
  % DDD
  % PPP

  Nrows = sim_in.Nsec*Rs;
  Nframes = floor((Nrows-1)/Ns);
  Nbits = Nframes * Nbitsperframe;    % number of payload data bits

  Nr = Nbits/Nc;                      % Number of data rows to get Nbits total

  if verbose
    printf("Nc.....: %d\n", Nc);
    printf("Ns.....: %d (step size for pilots, Ns-1 data symbols between pilots)\n", Ns);
    printf("Nr.....: %d\n", Nr);
    printf("Nbits..: %d\n", Nbits);
  end

  % double check if Nbits fit neatly into carriers

  assert(Nbits/Nc == floor(Nbits/Nc), "Nbits/Nc must be an integer");
 
  printf("Nframes: %d\n", Nframes);

  Nrp = Nr + Nframes + 1;  % number of rows once pilots inserted
                           % extra row of pilots at end
  printf("Nrp....: %d (number of rows including pilots)\n", Nrp);

  % simulate for each Eb/No point ------------------------------------

  for nn=1:length(EbNodB)
    rand('seed',1);
    randn('seed',1);

    EbNo = 10 .^ (EbNodB(nn)/10);
    variance = 1/(M*EbNo/2);

    % generate tx bits and insert pilot rows and border cols

    tx_bits = []; tx_bits_np = [];
    for f=1:Nframes
      tx_bits = [tx_bits; ones(1,Nc+2)];
      for r=1:Nrowsperframe
        arowofbits = rand(1,Nc) > 0.5;
        tx_bits = [tx_bits; [1 arowofbits 1]];
        tx_bits_np = [tx_bits_np; arowofbits];
      end
    end      
    tx_bits = [tx_bits; [1 ones(1,Nc) 1]]; % final row of pilots
    [nr nc] = size(tx_bits);
    assert(nr == Nrp);

    % map to BPSK symbols

    tx_sym = 2*tx_bits - 1;

#{
    % upsample to rate Fs using Zero Order Hold (ZOH)

    tx_sym_os = zeros(Nrp*M, Nc+2);
    for r=1:Nrp
      for rr=(r-1)*M+1:r*M
        tx_sym_os(rr,:) = tx_sym(r,:);
      end
    end
#}

    % OFDM up conversion and upsampling to rate Fs

    w = (0:Nc+1)*2*pi*Rs/Fs;
    W = zeros(M,Nc+2);
    for c=1:Nc+2
      W(:,c) = exp(j*w(c)*(0:M-1));
    end

    Nsam = Nrp*M;
    tx = zeros(Nrp*M,1);
    for r=1:Nrp
      for c=1:Nc+2
        acarrier = tx_sym(r,c) * W(:,c);
        tx((r-1)*M+1:r*M) += acarrier/M;
      end
    end
        
    rx = tx;

    noise = sqrt(variance)*(0.5*randn(Nrp*M,1) + j*0.5*randn(Nrp*M,1));
    rx += noise;

    % downconvert, downsample and integrate using OFDM

    rx_sym = zeros(Nrp, Nc+2);
    for r=1:Nrp
      for c=1:Nc+2
        acarrier = rx((r-1)*M+1:r*M) .* conj(W(:,c));
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

    % remove pilots to give us just data symbols
      
    rx_np = [];
    for r=1:Nrp
      if mod(r-1,Ns) != 0
        rx_np = [rx_np; rx_corr(r,2:Nc+1)];
      end
    end

    % calculate BER stats as a block, after pilots extracted

    rx_bits_np = real(rx_np) > 0;
    errors = xor(tx_bits_np, rx_bits_np);
    Nerrs = sum(sum(errors));

    printf("EbNodB: %3.2f BER: %4.3f Nbits: %d Nerrs: %d\n", EbNodB(nn), Nerrs/Nbits, Nbits, Nerrs);

    if verbose
      figure(1)
      plot(real(tx))
      figure(2)
      Tx = abs(fft(tx.*hanning(Nsam)));
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
        plot(abs(hf_model(:,2:Nc+1)));
      end

      figure(5); clf;
      plot(phase_est_log(:,2:Nc+1),'+', 'markersize', 10); 
      hold on; 
      plot(phase_est_pilot_log(:,2:Nc+1),'g+', 'markersize', 5); 
      if sim_in.hf_en
        plot(angle(hf_model(:,2:Nc+1)));
      end
      axis([1 Nrp -pi pi]);
    end

    sim_out.ber(nn) = sum(Nerrs)/Nbits; 
    sim_out.pilot_overhead = 10*log10(Ns/(Ns-1));
  end
endfunction


% Plot BER against Eb/No curves at various pilot insertion rates Ns
% using the HF multipath channel.  Second set of curves includes Eb/No
% loss for pilot insertion, so small Ns means better tracking of phase
% but large pilot insertion loss

% Target operating point Eb/No is 6dB, as this is where our rate 1/2
% LDPC code gives good results (10% PER, 1% BER).  However this means
% the Eb/No at the input is 10*log(1/2) or 3dB less, so we need to
% make sure phase est works at Eb/No = 6 - 3 = 3dB

function run_curves_hf
  sim_in.Nc = 8;
  sim_in.Ns = 8;
  sim_in.Nsec = 240;
  sim_in.EbNodB = 1:8;
  sim_in.verbose = 0;
  sim_in.hf_en = 0;

  EbNoLin = 10.^(sim_in.EbNodB/10);
  hf_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));

  hf_sim = run_sim(sim_in);

  figure(4); clf;
  semilogy(sim_in.EbNodB, hf_theory,'b+-;HF theory;');
  hold on;
  semilogy(sim_in.EbNodB, hf.ber,'c+-;Ns=8 HFsim;');
  hold off;
  axis([1 8 4E-2 2E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
end


function run_single
  sim_in.Nc = 8;
  sim_in.Ns = 8;
  sim_in.Nsec = 10;
  sim_in.EbNodB = 4;
  sim_in.verbose = 1;
  sim_in.hf_en = 0;

  run_sim(sim_in);
end


format;
more off;


run_single
%run_curves_hf




