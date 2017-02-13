% hf_modem_curves
% David Rowe Feb 2017
%
% Ideal implementations of a bunch of different HF modems, used to
% generate plots for a blog post.

#{
  [X] ideal AWGN/HF curves
  [X] exp AWGN QPSK curves
  [X] exp AWGN DQPSK curves
  [X] exp HF channel model
  [ ] diversity
  [ ] COHPSK frames
      + would require multiple carriers
      + filtering or OFDM
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


% Rate Rs modem simulation model -------------------------------------------------------

function sim_out = ber_test(sim_in)
    bps     = 2; % two bits/symbol for QPSK

    verbose = sim_in.verbose;
    EbNovec = sim_in.EbNovec;
    hf_en   = sim_in.hf_en;

    % user can supply number of bits per point to get good results
    % at high Eb/No

    if length(sim_in.nbits) > 1
      nbitsvec = sim_in.nbits;
      nbitsvec += 100 - mod(nbitsvec,100);  % round up to nearest 100
    else
      nbitsvec(1:length(EbNovec)) = sim_in.nbits;
    end

    % init HF model

    if hf_en
      % some typical values
      Rs = 50; dopplerSpreadHz = 1; path_delay = 0.001/Rs;

      nsymb = max(nbitsvec)/2;
      spread1 = doppler_spread(dopplerSpreadHz, Rs, nsymb);
      spread2 = doppler_spread(dopplerSpreadHz, Rs, nsymb);
      hf_gain = 1.0/sqrt(var(spread1)+var(spread2));
      %printf("nsymb: %d lspread1: %d\n", nsymb, length(spread1));
    end

    for ne = 1:length(EbNovec)

        % work out noise power -------------

        EbNodB = EbNovec(ne);
        EsNodB = EbNodB + 10*log10(bps);
        EsNo = 10^(EsNodB/10);
        variance = 1/EsNo;
        nbits = nbitsvec(ne);
        nsymb = nbits/bps;

        % modulator ------------------------

        tx_bits = rand(1,nbits) > 0.5;
        tx_symb = [];
        prev_tx_symb = 1;
        for s=1:nsymb
          atx_symb = qpsk_mod(tx_bits(2*s-1:2*s));
          if sim_in.dqpsk
            atx_symb *= prev_tx_symb;
            prev_tx_symb = atx_symb;
          end
          tx_symb = [tx_symb atx_symb];
        end

        % channel ---------------------------

        rx_symb = tx_symb;

        if hf_en

          % simplified rate Rs simulation model that doesn't include
          % ISI, just freq filtering.  We assume perfect phase estimation
          % so it's just amplitude distortion.
 
          for s=1:nsymb 
            hf_model = hf_gain*(spread1(s) + exp(-j*path_delay)*spread2(s));
            rx_symb(s) = rx_symb(s)*abs(hf_model);
          end
        end

        % variance is noise power, which is divided equally between real and
        % imag components of noise

        noise = sqrt(variance*0.5)*(randn(1,nsymb) + j*randn(1,nsymb));
        rx_symb += noise;

        % demodulator ------------------------------------------

        rx_bits = [];
        prev_rx_symb = 1;
        for s=1:nsymb
          arx_symb = rx_symb(s);
          if sim_in.dqpsk
            tmp = arx_symb;
            arx_symb *= prev_rx_symb';
            prev_rx_symb = tmp;
          end
          two_bits = qpsk_demod(arx_symb);
          rx_bits = [rx_bits two_bits];
        end
        
        % count errors -----------------------------------------

        error_pattern = xor(tx_bits, rx_bits);
        nerrors = sum(error_pattern);
        bervec(ne) = nerrors/nbits;
        if verbose
          printf("EbNodB: %3.1f nerrors: %5d ber: %4.3f\n", EbNodB, nerrors, bervec(ne));
          if verbose == 2
            figure(2); clf;
            plot(rx_symb,'+');
          end
        end
    end

    sim_out.bervec = bervec;
endfunction


% -------------------------------------------------------------


function run_single
    sim_in.verbose = 2;
    sim_in.nbits   = 10000;
    sim_in.EbNovec = 8;
    sim_in.dqpsk   = 0;
    sim_in.hf_en   = 1;

    sim_qpsk = ber_test(sim_in);
endfunction


function run_curves
    sim_in.verbose = 1;
    sim_in.EbNovec = 0:10;
    sim_in.dqpsk   = 0;
    sim_in.hf_en   = 0;

    % AWGN -----------------------------

    ber_awgn_theory = 0.5*erfc(sqrt(10.^(sim_in.EbNovec/10)));
    sim_in.nbits    = min(1E5, floor(500 ./ ber_awgn_theory));

    sim_qpsk = ber_test(sim_in);
    sim_in.dqpsk = 1;
    sim_dqpsk = ber_test(sim_in);
       
    % HF -----------------------------

    hf_sim_in = sim_in; hf_sim_in.dqpsk = 0; hf_sim_in.hf_en = 1;
    hf_sim_in.EbNovec = 0:16;

    EbNoLin = 10.^(hf_sim_in.EbNovec/10);
    ber_hf_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));

    hf_sim_in.nbits = min(1E5, floor(500 ./ ber_hf_theory));
    sim_qpsk_hf = ber_test(hf_sim_in);

    hf_sim_in.dqpsk = 1; 
    sim_dqpsk_hf = ber_test(hf_sim_in);

    % Plot results --------------------

    figure(1); clf;
    semilogy(sim_in.EbNovec, ber_awgn_theory,'r+-;QPSK AWGN theory;','markersize', 10, 'linewidth', 2)
    hold on;
    semilogy(hf_sim_in.EbNovec, ber_hf_theory,'r+-;QPSK HF theory;','markersize', 10, 'linewidth', 2)
    semilogy(sim_in.EbNovec, sim_qpsk.bervec,'g+-;QPSK AWGN simulated;','markersize', 10, 'linewidth', 2)
    semilogy(sim_in.EbNovec, sim_dqpsk.bervec,'b+-;DQPSK AWGN simulated;','markersize', 10, 'linewidth', 2)
    semilogy(hf_sim_in.EbNovec, sim_qpsk_hf.bervec,'c+-;QPSK HF simulated;','markersize', 10, 'linewidth', 2)
    semilogy(hf_sim_in.EbNovec, sim_dqpsk_hf.bervec,'k+-;DQPSK HF simulated;','markersize', 10, 'linewidth', 2)
    hold off;
    xlabel('Eb/N0')
    ylabel('BER')
    grid("minor")
    axis([min(hf_sim_in.EbNovec) max(hf_sim_in.EbNovec) 1E-3 1])

endfunction

% -------------------------------------------------------------

more off;
rand('seed',1); randn('seed', 1);
run_curves
%run_single
