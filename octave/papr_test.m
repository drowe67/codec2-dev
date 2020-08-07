% papr_test.m
%
#{

  TODO:
    [X] measure BER
    [ ] curves with different clipping
    [ ] companding
    [ ] heat map type scatter diagram
    [ ] tx diversity
    [ ] think about how to meausre/objective
        + can we express BER as a function of peak power?
        + or maybe showing RMS power being boosted
        + take into account clipping and higher rms power on BER
#}

1;

function symbol = qpsk_mod(two_bits)
    two_bits_decimal = sum(two_bits .* [2 1]); 
    switch(two_bits_decimal)
        case (0) symbol =  1;
        case (1) symbol =  j;
        case (2) symbol = -j;
        case (3) symbol = -1;
    endswitch
endfunction

function two_bits = qpsk_demod(symbol)
    bit0 = real(symbol*exp(j*pi/4)) < 0;
    bit1 = imag(symbol*exp(j*pi/4)) < 0;
    two_bits = [bit1 bit0];
endfunction

function [ber papr] = run_sim(Nsym,EbNodB,plot_en=0, clip=0)
    Nc   = 8;
    M    = 160;   % number of samples in each symbol
    bps  = 2;     % two bits per symbol for QPSK

    for e=1:length(EbNodB)
        % generate a 2D array of QPSK symbols

        Nphases = 2^bps;
        tx_phases = pi/2*floor((rand(Nsym,Nc)*Nphases));
        tx_sym = exp(j*tx_phases);

        w = 2*pi/M;
        tx = [];

        % generate OFDM signal

        for s=1:Nsym
          atx = zeros(1,M);
          for c=1:Nc
            atx += exp(j*(0:M-1)*c*w)*tx_sym(s,c);
          end
          tx = [tx atx];
        end
        Nsam = length(tx);

        if (clip)
            tx(find(abs(tx)>6)) = 6;
        end
        
        % AWGN channel

        EsNodB = EbNodB(e) + 10*log10(bps);
        variance = M/(10^(EsNodB/10));
        noise = sqrt(variance/2)*randn(1,Nsam) + j*sqrt(variance/2)*randn(1,Nsam);
        rx = tx + noise;

        % demodulate

        rx_sym = zeros(Nsym,Nc);
        for s=1:Nsym
          st = (s-1)*M+1; en = s*M;
          for c=1:Nc
            rx_sym(s,c) = sum(exp(-j*(0:M-1)*c*w) .* rx(st:en))/M;
          end
        end

        % count bit errors

        Tbits = Terrs = 0;
        for s=1:Nsym
          for c=1:Nc
            tx_bits = qpsk_demod(tx_sym(s,c));
            rx_bits = qpsk_demod(rx_sym(s,c));
            Tbits += bps;
            Terrs += sum(xor(tx_bits,rx_bits));
          end
        end

        if plot_en
          figure(1); clf;
          plot(abs(tx(1:10*M)))
          figure(2); clf; hist(abs(tx),25)
          figure(3); clf;
          rx_sym = reshape(rx_sym, Nsym*Nc, 1);
          h_data = [real(rx_sym)'; imag(rx_sym)']';
          h = log10(hist3(h_data,[50 50])+1);
          colormap("default"); imagesc(h);
          figure(4); clf; plot(real(rx_sym), imag(rx_sym), '+'); axis([-2 2 -2 2])
        end

        papr = 20*log10(max(abs(tx))/mean(abs(tx)));
        ber(e) = Terrs/Tbits;
        printf("EbNodB: %3.1f PAPR: %5.2f Tbits: %6d Terrs: %6d BER: %5.3f\n", EbNodB(e), papr, Tbits, Terrs, ber(e))
    end  
end

#run_sim(1000,10,1, clip_en=1)

Nsym=1000;
EbNodB=3:8;
[ber1 papr1] = run_sim(Nsym,EbNodB);
[ber2 papr2] = run_sim(Nsym,EbNodB,0,clip_en=1);

figure(1); clf;
semilogy(EbNodB, ber1,'b'); hold on;
semilogy(EbNodB, ber2,'r'); hold off;
axis([min(EbNodB) max(EbNodB) 1E-4 1])
