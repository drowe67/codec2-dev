% hf_modem_curves
% David Rowe Feb 2017
%
% Ideal implementations of a bunch of different HF modems, used to
% generate plots for a blog post.

#{
  [ ] ideal AWGN/HF curves
  [ ] exp AWGN QPSK curves
  [ ] exp AWGN DQPSK curves
  [ ] exp HF channel model
  [ ] with COHPSK frames
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


function sim_out = ber_test(sim_in)
    bps     = 2; % two bits/symbol for QPSK

    verbose = sim_in.verbose;
    nbits   = sim_in.nbits;
    nsymb   = nbits/bps;
    EbNovec = sim_in.EbNovec;

    for ne = 1:length(EbNovec)
        EbNodB = EbNovec(ne);
        EsNodB = EbNodB + 10*log10(bps);
        EsNo = 10^(EsNodB/10);
    
        variance = 1/EsNo;

        tx_bits = rand(1,nbits) > 0.5;
        tx_symb = [];
        for s=1:nsymb
          atx_symb = qpsk_mod(tx_bits(2*s-1:2*s));
          tx_symb = [tx_symb atx_symb];
        end

        % variance is noise power, which is divided equally between real and
        % imag components of noise

        noise = sqrt(variance*0.5)*(randn(1,nsymb) + j*randn(1,nsymb));
        rx_symb = tx_symb + noise;

        rx_bits = [];
        for s=1:nsymb
          two_bits = qpsk_demod(rx_symb(s));
          rx_bits = [rx_bits two_bits];
        end
        
        error_pattern = xor(tx_bits, rx_bits);
        nerrors = sum(error_pattern);
        bervec(ne) = nerrors/nbits;
        if verbose
          printf("EbNodB: %3.1f nerrors: %5d ber: %4.3f\n", EbNodB, nerrors, bervec(ne));
        end
    end

    sim_out.bervec = bervec;
endfunction

% -------------------------------------------------------------

function run_curves
    sim_in.verbose = 1;
    sim_in.nbits   = 100000;
    sim_in.EbNovec = 0:7;

    ber_theory = 0.5*erfc(sqrt(10.^(sim_in.EbNovec/10)));
    sim_qpsk = ber_test(sim_in);

    figure(1); clf;
    semilogy(sim_in.EbNovec, ber_theory,'r+-;QPSK AWGN theory;','markersize', 10, 'linewidth', 2)
    hold on;
    semilogy(sim_in.EbNovec, sim_qpsk.bervec,'g+-;QPSK AWGN simulated;','markersize', 10, 'linewidth', 2)
    hold off;
    xlabel('Eb/N0')
    ylabel('BER')
    grid("minor")
    axis([min(sim_in.EbNovec) max(sim_in.EbNovec) 1E-3 1])

endfunction

% -------------------------------------------------------------

more off;
rand('seed',1); randn('seed', 1);
run_curves
