% papr_test.m
%
#{

  TODO:
    [X] measure BER
    [X] heat map type scatter diagram
    [ ] curves with different clipping thresholds
    [X] clipping
    [X] companding
    [ ] set threshold at CDF
    [ ] filter
        + complex or real?
        + might introduce phase shift ... make demodulation tricky
    [ ] way to put plots on top of each other
    [ ] tx diversity
    [ ] think about how to meausre/objective
        + can we express BER as a function of peak power?
        + or maybe showing RMS power being boosted
        + take into account clipping and higher rms power on BER
        + do we need to recalculate tx power?
    [ ] down the track ... other types of channels
        + does it mess up HF multipath much
    [ ] Tx a companded waveform and see if it really helps average power
    [ ] results
        + we can trade of PAPR with bandwidth, lower PAPR, more bandwidth
        + filtering bring sthe PAPR up again
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

function [ber papr_log] = run_sim(Nsym, EbNodB, plot_en=0, filt_en=0, method="")
    Nc   = 16;
    M    = 160;   % number of samples in each symbol
    bps  = 2;     % two bits per symbol for QPSK

    papr_log = [];
    for e=1:length(EbNodB)
        % generate a 2D array of QPSK symbols

        Nphases = 2^bps;
        tx_phases = pi/2*floor((rand(Nsym,Nc)*Nphases));
        tx_sym = exp(j*tx_phases);

        w = 2*pi/M; woff = (-Nc/2+0.5)*2*pi/M;
        tx = [];

        % generate OFDM signal

        for s=1:Nsym
          atx = zeros(1,M);
          for c=1:Nc
            atx += exp(j*(0:M-1)*(c*w+woff))*tx_sym(s,c);
          end
          tx = [tx atx];
        end
        Nsam = length(tx);

        tx_ = tx;
        if strcmp(method, "clip")
            ind = find(abs(tx) > 6);
            tx_(ind) = 6*exp(j*angle(tx(ind)));
        end
        if strcmp(method, "compand")
          tx_mag = interp1([0 1 3 5 6 20], [0 1 4 5.5 6 7], abs(tx), "pchip");
          tx_ = tx_mag.*exp(j*angle(tx));
        end

        if filt_en
          Nfilt=40;
          b = fir1(Nfilt,1.5*Nc/M);
          tx_ = filter(b,1,[tx_ zeros(1,Nfilt/2)]);
          tx_ = [tx_(Nfilt/2+1:end)];
        end
        
        % AWGN channel

        EsNodB = EbNodB(e) + 10*log10(bps);
        variance = M/(10^(EsNodB/10));
        noise = sqrt(variance/2)*randn(1,Nsam) + j*sqrt(variance/2)*randn(1,Nsam);
        rx = tx_ + noise;

        % demodulate

        rx_sym = zeros(Nsym,Nc);
        for s=1:Nsym
          st = (s-1)*M+1; en = s*M;
          for c=1:Nc
            rx_sym(s,c) = sum(exp(-j*(0:M-1)*(c*w+woff)) .* rx(st:en))/M;
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
          plot(abs(tx(1:5*M))); hold on; plot(abs(tx_(1:5*M))); hold off;
          axis([0 5*M 0 max(abs(tx))]);
          figure(2); clf; hist(abs(tx),25);
          % heat map type scatter plot
          figure(3); clf;
          rx_sym = reshape(rx_sym, Nsym*Nc, 1);
          h_data = [real(rx_sym)'; imag(rx_sym)']';
          h = log10(hist3(h_data,[50 50])+1);
          colormap("default"); imagesc(h);
          figure(4); clf; plot(real(rx_sym), imag(rx_sym), '+'); axis([-2 2 -2 2])
          figure(5); clf; Tx_ = 10*log10(abs(fft(tx_))); plot(fftshift(Tx_));
          mx = 10*ceil(max(Tx_)/10); axis([1 length(Tx_) mx-60 mx]);
        end

        papr1 = 20*log10(max(abs(tx))/mean(abs(tx)));
        papr2 = 20*log10(max(abs(tx_))/mean(abs(tx_)));
        papr_log = [papr_log papr2];
        ber(e) = Terrs/Tbits;
        printf("EbNodB: %3.1f PAPR: %5.2f %5.2f Tbits: %6d Terrs: %6d BER: %5.3f\n", EbNodB(e), papr1, papr2, Tbits, Terrs, ber(e))
    end

    papr = mean(papr_log);
end

run_sim(1000,100,plot_en=1, filt_en=1, "compand")

Nsym=1000;
EbNodB=3:8;
[ber1 papr1] = run_sim(Nsym,EbNodB, 0, filt_en=1);
[ber2 papr2] = run_sim(Nsym,EbNodB, 0, filt_en=1, "clip");
[ber3 papr3] = run_sim(Nsym,EbNodB, 0, filt_en=1, "compand");

figure(6); clf;
semilogy(EbNodB, ber1,sprintf('b;vanilla OFDM %3.1f;',papr1)); hold on;
semilogy(EbNodB, ber2,sprintf('r;clip %3.1f;',papr2)); 
semilogy(EbNodB, ber3,sprintf('g;compand %3.1f;',papr3)); hold off;
axis([min(EbNodB) max(EbNodB) 1E-4 1E-1]); grid;
