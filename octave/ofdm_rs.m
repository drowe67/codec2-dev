% ofdm_rs.m
% David Rowe Mar 2017
%
% Rate Rs BPSK/QPSK simulation to explore techniques for
% phase estimation over multiple carriers in HF channel.  Rate
% Rs version of ofdm_fs.m

#{
 TODO:
   [X] sim pilot based phase est using known symbols
   [X] test AWGN BER with averaging pilots from adj carriers
   [X] refactor to insert pilot rows
   [X] add border cols, not used for data
   [X] centre est on current carrier, extend to > 3
   [X] test single points
       + 1dB IL @ 6dB HF, 0.4 dB @ 2dB AWGN
   [X] try linear interpolation
   [X] try longer time windows
   [X] try combining mod stripping phase est inside frame
   [X] curves taking into account pilot losses
   [ ] remove border carriers, interpolate edge carrier
   [X] modify for QPSK
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
  Rs = 100;

  bps = sim_in.bps;
  EbNodB = sim_in.EbNodB;
  verbose = sim_in.verbose;
  hf_en = sim_in.hf_en;
  hf_phase = sim_in.hf_phase;
  phase_offset = sim_in.phase_offset;

  Ns = sim_in.Ns;          % step size for pilots
  Nc = sim_in.Nc;          % Number of cols, aka number of carriers

  Nbitsperframe = (Ns-1)*Nc*bps;
  Nrowsperframe = Nbitsperframe/(Nc*bps);
  if verbose
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

  Nrows = sim_in.Nsec*Rs;
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
  
  % set up HF model

  if hf_en

    % some typical values, or replace with user supplied

    dopplerSpreadHz = 1.0; path_delay = 1E-3*Rs;

    if isfield(sim_in, "dopplerSpreadHz") 
      dopplerSpreadHz = sim_in.dopplerSpreadHz;
    end
    if isfield(sim_in, "path_delay") 
      path_delay = sim_in.path_delay;
    end
    printf("Doppler Spread: %3.2f Hz Path Delay: %3.2f symbols\n", dopplerSpreadHz, path_delay);
    randn('seed',1);
    spread1 = doppler_spread(dopplerSpreadHz, Rs, sim_in.Nsec*Rs*1.1);
    spread2 = doppler_spread(dopplerSpreadHz, Rs, sim_in.Nsec*Rs*1.1);

    % sometimes doppler_spread() doesn't return exactly the number of samples we need
 
    assert(length(spread1) >= Nrp, "not enough doppler spreading samples");
    assert(length(spread2) >= Nrp, "not enough doppler spreading samples");
  end
  
  % construct an artificial phase countour for testing, linear across freq and time

  if sim_in.phase_test
    phase_test = ones(Nrp, Nc+2);
    for r=1:Nrp
      for c=1:Nc+2
        phase_test(r,c) = -pi/2 + c*pi/(Nc+2) + r*0.01*2*pi;
        phase_test(r,c) = phase_test(r,c) - 2*pi*floor((phase_test(r,c)+pi)/(2*pi));
      end
    end
  end

  % simulate for each Eb/No point ------------------------------------

  for nn=1:length(EbNodB)
    rand('seed',1);
    randn('seed',1);

    EsNo = bps * (10 .^ (EbNodB(nn)/10));
    variance = 1/(EsNo/2);
    noise = sqrt(variance)*(0.5*randn(Nrp,Nc+2) + j*0.5*randn(Nrp,Nc+2));

    % generate tx bits

    tx_bits = rand(1,Nbits) > 0.5;

    % map to symbols 

    if bps == 1
      tx_symb = 2*tx_bits - 1;
    end
    if bps == 2
      for s=1:Nbits/bps
        tx_symb(s) = qpsk_mod(tx_bits(2*(s-1)+1:2*s));
      end
    end

    % place symbols in multi-carrier frame with pilots and boundary carriers

    tx = []; s = 1;
    for f=1:Nframes
      aframe = zeros(Nrowsperframe,Nc+2);
      aframe(1,:) = 1;
      for r=1:Nrowsperframe
        arowofsymbols = tx_symb(s:s+Nc-1);
        s += Nc;
        aframe(r+1,2:Nc+1) = arowofsymbols;
      end
      tx = [tx; aframe];
    end      
    tx = [tx; ones(1,Nc+2)];  % final row of pilots
    [nr nc] = size(tx);
    assert(nr == Nrp);
    
    rx = tx * exp(j*phase_offset);

    if sim_in.phase_test
      rx = rx .* exp(j*phase_test);
    end

    if hf_en

      % simplified rate Rs simulation model that doesn't include
      % ISI, just freq filtering.
      
      % Note Rs carrier spacing, sample rate is Rs

      hf_model = zeros(Nr,Nc+2); phase_est = zeros(Nr,Nc);
      for r=1:Nrp
        for c=1:Nc+2
          w = 2*pi*c*Rs/Rs; 
          hf_model(r,c) = spread1(r) + exp(-j*w*path_delay)*spread2(r);
        end
        
        if hf_phase
          rx(r,:) = rx(r,:) .* hf_model(r,:);
        else
          rx(r,:) = rx(r,:) .* abs(hf_model(r,:));
        end
      end

      % normalise power over HF simulation run

      p = sum(var(rx));
      rx *= sqrt(Nc/p);
    end

    rx += noise;

    % pilot based phase est, we use known tx symbols as pilots ----------

    rx_corr = rx;

    if sim_in.pilot_phase_est

      % est phase from pilots either side of data symbols
      % adjust phase of data symbol
      % demodulate and count errors of just data
 
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
          aphase_est_pilot_rect1 = sum(rx(r,cr)*tx(r,cr)');
          aphase_est_pilot_rect2 = sum(rx(r+Ns,cr)*tx(r+Ns,cr)');

          % optionally use next step of pilots in past and future

          if sim_in.pilot_wide
            if r > Ns+1
              aphase_est_pilot_rect1 += sum(rx(r-Ns,cr)*tx(r-Ns,cr)');
            end
            if r < Nrp - 2*Ns
              aphase_est_pilot_rect2 += sum(rx(r+2*Ns,cr)*tx(r+2*Ns,cr)');
            end
          end

          % correct phase offset using phase estimate

          for rr=r+1:r+Ns-1
            a = b = 1;
            if sim_in.pilot_interp
              b = (rr-r)/Ns; a = 1 - b;
            end
            %printf("rr: %d a: %4.3f b: %4.3f\n", rr, a, b);
            aphase_est_pilot = angle(a*aphase_est_pilot_rect1 + b*aphase_est_pilot_rect2);
            phase_est_pilot_log(rr,c) = aphase_est_pilot;
            rx_corr(rr,c) = rx(rr,c) * exp(-j*aphase_est_pilot);
          end

          if sim_in.stripped_phase_est
            % Optional modulation stripping feed fwd phase estimation, to refine
            % pilot-based phase estimate.  Doing it after pilot based phase estimation
            % means we don't need to deal with ambiguity, which is difficult to handle
            % in low SNR channels.

            % Use vector of 7 symbols around current data symbol.  We could use a 2D
            % window if we can work out how best to correct with pilot-est and avoid
            % ambiguities

            for rr=r+1:r+Ns-1

              % extract a matrix of nearby samples with pilot-based offset removed
 
              amatrix = rx(max(1,rr-3):min(Nrp,rr+3),c) .* exp(-j*aphase_est_pilot);

              % modulation strip and est phase

              stripped = abs(amatrix) .* exp(j*2*angle(amatrix));
              aphase_est_stripped = angle(sum(sum(stripped)))/2;
              phase_est_stripped_log(rr,c) = aphase_est_stripped;

              % correct rx symbols based on both phase ests

              phase_est_log(rr,c) = angle(exp(j*(aphase_est_pilot+aphase_est_stripped)));
              rx_corr(rr,c) = rx(rr,c) * exp(-j*phase_est_log(rr,c));
            end       
          end % sim_in.stripped_phase_est

        end % r=1:Ns:Nrp-Ns

      end % c=2:Nc+1
    end % sim_in.pilot_phase_est


    if isfield(sim_in, "ml_pd") && sim_in.ml_pd

      % Bill's ML with pilots phase detector, does phase est and demodulation 

      rx_bits = []; rx_np = [];
      aframeofbits = zeros(Ns-1, Nc);
      for r=1:Ns:Nrp-Ns

        % demodulate this frame, ML operates carrier by carrier 

        for c=2:Nc+1
          arxcol = rx(r:r+Ns, c);
          arxcol(1) = rx(r, c-1) + rx(r, c+1);
          arxcol(Ns+1) = rx(r+Ns, c-1) + rx(r+Ns, c+1);
          [acolofbits aphase_est] = ml_pd(rot90(arxcol),  bps, [1 Ns+1]);
          aframeofbits(:,c-1) = xor(acolofbits, ones(1,Ns-1));
          rx_np = [rx_np rot90(arxcol) .* exp(-j*aphase_est)];
        end
        
        % unpack from frame into linear array of bits

        for rr=1:Ns-1
          rx_bits = [rx_bits aframeofbits(rr,:)];
        end
 
      end
    else

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
    end
    %tx_bits
    %rx_bits
    assert(length(rx_bits) == Nbits);

    %phase_test
    %phase_est_log

    % calculate BER stats as a block, after pilots extracted

    errors = xor(tx_bits, rx_bits);
    Nerrs = sum(errors);

    printf("EbNodB: %3.2f BER: %5.4f Nbits: %d Nerrs: %d\n", EbNodB(nn), Nerrs/Nbits, Nbits, Nerrs);

    if verbose
      figure(1); clf; 
      plot(rx_np,'+');
      axis([-2 2 -2 2]);

      if hf_en
        figure(2); clf; 
        plot(abs(hf_model(:,2:Nc+1)));
      end

      if sim_in.pilot_phase_est
        figure(3); clf;
        plot(phase_est_log(:,2:Nc+1),'+', 'markersize', 10); 
        hold on; 
        plot(phase_est_pilot_log(:,2:Nc+1),'g+', 'markersize', 5); 
        if sim_in.stripped_phase_est
          plot(phase_est_stripped_log(:,2:Nc+1),'ro', 'markersize', 5); 
        end
        if sim_in.hf_en && sim_in.hf_phase
          plot(angle(hf_model(:,2:Nc+1)));
        end
        if sim_in.phase_test
          plot(phase_test(:,2:Nc+1));
        end
        axis([1 Nrp -pi pi]);
      end
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
  sim_in.Nc = 7;
  sim_in.Ns = 5;
  sim_in.Nsec = 240;
  sim_in.EbNodB = 1:8;
  sim_in.verbose = 0;
  sim_in.pilot_phase_est = 0;
  sim_in.pilot_wide = 1;
  sim_in.pilot_interp = 0;
  sim_in.stripped_phase_est = 0;
  sim_in.phase_offset = 0;
  sim_in.phase_test = 0;
  sim_in.hf_en = 1;
  sim_in.hf_phase = 0;

  sim_in.Ns = 5;
  hf_ref_Ns_5_no_phase = run_sim(sim_in);
  sim_in.Ns = 9;
  hf_ref_Ns_9_no_phase = run_sim(sim_in);

  sim_in.hf_phase = 1;
  sim_in.pilot_phase_est = 1;

  sim_in.Ns = 5;
  hf_Ns_5 = run_sim(sim_in);

  sim_in.Ns = 9;
  hf_Ns_9 = run_sim(sim_in);

  sim_in.Ns = 17;
  hf_Ns_17 = run_sim(sim_in);

  figure(4); clf;
  semilogy(sim_in.EbNodB, hf_ref_Ns_5_no_phase.ber,'b+-;Ns=5 HF ref no phase;');
  hold on;
  semilogy(sim_in.EbNodB, hf_ref_Ns_9_no_phase.ber,'c+-;Ns=9 HF ref no phase;');
  semilogy(sim_in.EbNodB, hf_Ns_5.ber,'g+--;Ns=5;');
  semilogy(sim_in.EbNodB + hf_Ns_5.pilot_overhead, hf_Ns_5.ber,'go-;Ns=5 with pilot overhead;');
  semilogy(sim_in.EbNodB, hf_Ns_9.ber,'r+--;Ns=9;');
  semilogy(sim_in.EbNodB + hf_Ns_9.pilot_overhead, hf_Ns_9.ber,'ro-;Ns=9 with pilot overhead;');
  semilogy(sim_in.EbNodB, hf_Ns_17.ber,'k+--;Ns=17;');
  semilogy(sim_in.EbNodB + hf_Ns_17.pilot_overhead, hf_Ns_17.ber,'ko-;Ns=17 with pilot overhead;');
  hold off;
  axis([1 8 4E-2 2E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
  title('HF Multipath 1Hz Doppler 1ms delay');

end


% Generate HF curves for some alternative, experimental methods tested
% during development, such as interpolation, refinements using
% modulation stripping, narrow window.

function run_curves_hf_alt
  sim_in.Nc = 7;
  sim_in.Ns = 5;
  sim_in.Nsec = 60;
  sim_in.EbNodB = 1:8;
  sim_in.verbose = 0;
  sim_in.pilot_phase_est = 0;
  sim_in.pilot_wide = 1;
  sim_in.pilot_interp = 0;
  sim_in.stripped_phase_est = 0;
  sim_in.phase_offset = 0;
  sim_in.phase_test = 0;
  sim_in.hf_en = 1;
  sim_in.hf_phase = 0;

  sim_in.Ns = 9;
  hf_ref_Ns_9_no_phase = run_sim(sim_in);

  sim_in.hf_phase = 1;
  sim_in.pilot_phase_est = 1;
  hf_Ns_9 = run_sim(sim_in);

  sim_in.stripped_phase_est = 1;
  hf_Ns_9_stripped = run_sim(sim_in);

  sim_in.stripped_phase_est = 0;
  sim_in.pilot_wide = 0;
  hf_Ns_9_narrow = run_sim(sim_in);

  sim_in.pilot_wide = 1;
  sim_in.pilot_interp = 1;
  hf_Ns_9_interp = run_sim(sim_in);

  figure(6); clf;
  semilogy(sim_in.EbNodB, hf_ref_Ns_9_no_phase.ber,'c+-;Ns=9 HF ref no phase;');
  hold on;
  semilogy(sim_in.EbNodB, hf_Ns_9.ber,'r+--;Ns=9;');
  semilogy(sim_in.EbNodB, hf_Ns_9_stripped.ber,'g+--;Ns=9 stripped refinement;');
  semilogy(sim_in.EbNodB, hf_Ns_9_narrow.ber,'b+--;Ns=9 narrow;');
  semilogy(sim_in.EbNodB, hf_Ns_9_interp.ber,'k+--;Ns=9 interp;');
  hold off;
  axis([1 8 4E-2 2E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
  title('HF Multipath 1Hz Doppler 1ms delay');

end


% Generate HF curves for fixed Ns but different HF channels.

function run_curves_hf_channels
  sim_in.Nc = 7;
  sim_in.Ns = 9;
  sim_in.Nsec = 240;
  sim_in.EbNodB = 1:8;
  sim_in.verbose = 0;
  sim_in.pilot_phase_est = 0;
  sim_in.pilot_wide = 1;
  sim_in.pilot_interp = 0;
  sim_in.stripped_phase_est = 0;
  sim_in.phase_offset = 0;
  sim_in.phase_test = 0;
  sim_in.hf_en = 1;
  sim_in.hf_phase = 0;

  hf_Ns_9_1hz_1ms_no_phase = run_sim(sim_in);

  sim_in.hf_phase = 1;
  sim_in.pilot_phase_est = 1;
  hf_Ns_9_1hz_1ms = run_sim(sim_in);

  Rs = 100;

  sim_in.dopplerSpreadHz = 1.0; 
  sim_in.path_delay = 500E-6*Rs;
  hf_Ns_9_1hz_500us = run_sim(sim_in);

  sim_in.dopplerSpreadHz = 1.0; 
  sim_in.path_delay = 2E-3*Rs;
  hf_Ns_9_1hz_2ms = run_sim(sim_in);

  sim_in.dopplerSpreadHz = 2.0; 
  sim_in.path_delay = 1E-3*Rs;
  hf_Ns_9_2hz_1ms = run_sim(sim_in);

  sim_in.dopplerSpreadHz = 2.0; 
  sim_in.path_delay = 1E-3*Rs;
  hf_Ns_9_2hz_2ms = run_sim(sim_in);

  sim_in.dopplerSpreadHz = 2.0; 
  sim_in.path_delay = 2E-3*Rs;
  hf_Ns_9_2hz_2ms = run_sim(sim_in);

  sim_in.dopplerSpreadHz = 4.0; 
  sim_in.path_delay = 1E-3*Rs;
  hf_Ns_9_4hz_1ms = run_sim(sim_in);

  figure(6); clf;
  semilogy(sim_in.EbNodB, hf_Ns_9_1hz_1ms_no_phase.ber,'c+-;Ns=9 1Hz 1ms ref no phase;');
  hold on;
  semilogy(sim_in.EbNodB, hf_Ns_9_1hz_500us.ber,'k+-;Ns=9 1Hz 500us;');
  semilogy(sim_in.EbNodB, hf_Ns_9_1hz_1ms.ber,'r+-;Ns=9 1Hz 1ms;');
  semilogy(sim_in.EbNodB, hf_Ns_9_1hz_2ms.ber,'bo-;Ns=9 1Hz 2ms;');
  semilogy(sim_in.EbNodB, hf_Ns_9_2hz_1ms.ber,'g+-;Ns=9 2Hz 1ms;');
  semilogy(sim_in.EbNodB, hf_Ns_9_2hz_2ms.ber,'mo-;Ns=9 2Hz 2ms;');
  semilogy(sim_in.EbNodB, hf_Ns_9_4hz_1ms.ber,'c+-;Ns=9 4Hz 1ms;');
  hold off;
  axis([1 8 4E-2 2E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
  title('HF Multipath Ns = 9');

end


% AWGN curves for BPSK and QPSK.  Coded Eb/No operating point is 2dB,
% so raw BER for rate 1/2 will be -1dB

function run_curves_awgn_bpsk_qpsk
  sim_in.Nc = 7;
  sim_in.Ns = 7;
  sim_in.Nsec = 30;
  sim_in.verbose = 0;
  sim_in.pilot_phase_est = 0;
  sim_in.pilot_wide = 1;
  sim_in.pilot_interp = 0;
  sim_in.stripped_phase_est = 1;
  sim_in.phase_offset = 0;
  sim_in.phase_test = 0;
  sim_in.hf_en = 0;
  sim_in.hf_phase = 0;

  sim_in.EbNodB = -3:5;

  ber_awgn_theory = 0.5*erfc(sqrt(10.^(sim_in.EbNodB/10)));

  sim_in.bps = 1;
  awgn_bpsk = run_sim(sim_in);
  sim_in.bps = 2;
  awgn_qpsk = run_sim(sim_in);

  figure(5); clf;
  semilogy(sim_in.EbNodB, ber_awgn_theory,'b+-;AWGN Theory;');
  hold on;
  semilogy(sim_in.EbNodB, awgn_bpsk.ber,'g+-;Ns=7 BPSK;');
  semilogy(sim_in.EbNodB + awgn_bpsk.pilot_overhead, awgn_bpsk.ber,'go-;Ns=7 BPSK with pilot overhead;');
  semilogy(sim_in.EbNodB, awgn_qpsk.ber,'r+-;Ns=7 QPSK;');
  semilogy(sim_in.EbNodB + awgn_qpsk.pilot_overhead, awgn_qpsk.ber,'ro-;Ns=7 QPSK with pilot overhead;');
  hold off;
  axis([-3 5 4E-3 2E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
  title('AWGN');
end


% HF multipath curves for BPSK and QPSK.  Coded operating point is about 3dB

function run_curves_hf_bpsk_qpsk
  sim_in.Nc = 7;
  sim_in.Ns = 7;
  sim_in.Nsec = 120;
  sim_in.verbose = 0;
  sim_in.pilot_phase_est = 1;
  sim_in.pilot_wide = 1;
  sim_in.pilot_interp = 0;
  sim_in.stripped_phase_est = 0;
  sim_in.phase_offset = 0;
  sim_in.phase_test = 0;
  sim_in.hf_en = 1;
  sim_in.hf_phase = 1;

  sim_in.EbNodB = 1:8;

  EbNoLin = 10.^(sim_in.EbNodB/10);
  hf_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));

  sim_in.bps = 1;
  hf_bpsk = run_sim(sim_in);
  sim_in.bps = 2;
  hf_qpsk = run_sim(sim_in);

  figure(5); clf;
  semilogy(sim_in.EbNodB, hf_theory,'b+-;HF Theory;');
  hold on;
  semilogy(sim_in.EbNodB, hf_bpsk.ber,'g+-;Ns=7 BPSK;');
  semilogy(sim_in.EbNodB + hf_bpsk.pilot_overhead, hf_bpsk.ber,'go-;Ns=7 BPSK with pilot overhead;');
  semilogy(sim_in.EbNodB, hf_qpsk.ber,'r+-;Ns=7 QPSK;');
  semilogy(sim_in.EbNodB + hf_qpsk.pilot_overhead, hf_qpsk.ber,'ro-;Ns=7 QPSK with pilot overhead;');
  hold off;
  axis([1 8 4E-3 2E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
  title('HF Multipath');
end


% AWGN curves for BPSK using 3 carrier 2D matrix pilot and ML pilot

function run_curves_awgn_ml
  sim_in.bps = 1;
  sim_in.Nc = 7;
  sim_in.Ns = 7;
  sim_in.Nsec = 10;
  sim_in.verbose = 0;
  sim_in.pilot_phase_est = 1;
  sim_in.pilot_wide = 1;
  sim_in.pilot_interp = 0;
  sim_in.stripped_phase_est = 1;
  sim_in.phase_offset = 0;
  sim_in.phase_test = 0;
  sim_in.hf_en = 0;
  sim_in.hf_phase = 0;

  sim_in.EbNodB = -3:5;

  ber_awgn_theory = 0.5*erfc(sqrt(10.^(sim_in.EbNodB/10)));

  awgn_2d = run_sim(sim_in);
  sim_in.pilot_phase_est = 0;
  sim_in.ml_pd = 1;
  awgn_ml = run_sim(sim_in);

  figure(5); clf;
  semilogy(sim_in.EbNodB, ber_awgn_theory,'b+-;AWGN Theory;');
  hold on;
  semilogy(sim_in.EbNodB, awgn_2d.ber,'g+-;Ns=7 3 carrier pilot BPSK;');
  semilogy(sim_in.EbNodB, awgn_ml.ber,'ro-;Ns=7 ML pilot BPSK;');
  hold off;
  axis([-3 5 4E-3 5E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
  title('AWGN');
end


% HF multipath curves for ML

function run_curves_hf_ml
  sim_in.bps = 1;
  sim_in.Nc = 7;
  sim_in.Ns = 14;
  sim_in.Nsec = 120;
  sim_in.verbose = 0;
  sim_in.pilot_phase_est = 1;
  sim_in.pilot_wide = 1;
  sim_in.pilot_interp = 0;
  sim_in.stripped_phase_est = 0;
  sim_in.phase_offset = 0;
  sim_in.phase_test = 0;
  sim_in.hf_en = 1;
  sim_in.hf_phase = 1;

  sim_in.EbNodB = 1:8;

  EbNoLin = 10.^(sim_in.EbNodB/10);
  hf_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));

  hf_2d = run_sim(sim_in);
  sim_in.pilot_phase_est = 0;
  sim_in.ml_pd = 1;
  hf_ml = run_sim(sim_in);

  figure(7); clf;
  semilogy(sim_in.EbNodB, hf_theory,'b+-;HF Theory;');
  hold on;
  semilogy(sim_in.EbNodB, hf_2d.ber,'g+-;Ns=7 3 carrier pilot BPSK;');
  semilogy(sim_in.EbNodB, hf_ml.ber,'ro-;Ns=7 ML pilot BPSK;');
  hold off;
  axis([1 8 4E-3 2E-1])
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid; grid minor on;
  legend('boxoff');
  title('HF Multipath');
end



function run_single
  sim_in.bps = 2;
  sim_in.Nsec = 60;
  sim_in.Nc = 16;
  sim_in.Ns = 8;
  sim_in.EbNodB = 6;
  sim_in.verbose = 1;

  sim_in.pilot_phase_est = 1;
  sim_in.pilot_wide = 1;
  sim_in.pilot_interp = 0;
  sim_in.stripped_phase_est = 0;
  sim_in.ml_pd = 0;

  sim_in.phase_offset = 0;
  sim_in.phase_test = 0;
  sim_in.hf_en = 1;
  sim_in.hf_phase = 1;
  sim_in.path_delay = 0;

  run_sim(sim_in);

  EbNoLin = 10.^(sim_in.EbNodB/10);
  hf_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));
  printf("HF theory: %5.4f\n", hf_theory);
end


format;
more off;

run_single
%run_curves_hf_bpsk_qpsk
%run_curves_hf_channels
%run_curves_hf_ml




