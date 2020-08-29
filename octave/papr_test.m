% papr_test.m
%
#{

  TODO:
    [X] measure BER
    [X] heat map type scatter diagram
    [X] clipping
    [X] companding
    [X] curves with different clipping thresholds
    [X] filter
        + complex or real?
        + might introduce phase shift ... make demodulation tricky
    [X] way to put plots on top of each other
    [X] tx diversity
    [ ] PAPR versus carriers for random data
    [ ] set threshold at CDF %
    [ ] plot compander
    [ ] think about how to measure/objective
        + can we express BER as a function of peak power?
        + or maybe showing RMS power being boosted
        + take into account clipping and higher rms power on BER
        + do we need to recalculate tx power?
    [ ] down the track ... other types of channels
        + does it mess up HF multipath much
        + impact likely less
    [ ] Tx a companded waveform and see if it really helps average power
    [ ] results/discussion
        + non-linear technqiues spread the energy in freq
        + we can trade of PAPR with bandwidth, lower PAPR, more bandwidth
        + filtering brings the PAPR up again
        + PAPR helps SNR, but multipath is a much tougher problem
        + non-linear technqiques will mess with QAM more (I think - simulate?)
        + so we may hit a wall at high data rates
        + are people already doing this by default, e.g. some slight compression in radio?  OTA test rqd,
          some sort of controlled test, like clipping rate
        + diversity idea doesn't seem to work in term of PAPR - but - 1 2Nc waveform will be much more
          robust to multipath (and multipath first), so any PAPR reduction will go further.
        + bigger improvements at higher Nc, e.g. 6dB at Nc=17.  Will need to test on real PAs
    [ ] what will happen when LDPC demodulated?
        + will decoder struggle due to non-linearity?
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

% "Genie" OFDM modem simulation that assumes ideal sync

function [ber papr] = run_sim(Nsym, EbNodB, channel='awgn', plot_en=0, filt_en=0, method="", threshold=6)
    Nc   = 8;
    M    = 160;    % number of samples in each symbol
    bps  = 2;      % two bits per symbol for QPSK
    Ncp  = 16;     % cyclic prefix samples
    Fs   = 8000;
    
    phase_est = 1; % perform phase estimation/correction
    timing = Ncp;

    if strcmp(method,"diversity")
      Nd = 2; gain = 1/sqrt(2);
    else
      Nd = 1; gain = 1.0;
    end
    
    if strcmp(channel,'multipath')
      dopplerSpreadHz = 1; path_delay = Ncp/2;
      Nsam = floor(Nsym*(M+Ncp)*1.1);
      spread1 = doppler_spread(dopplerSpreadHz, Fs, Nsam);
      spread2 = doppler_spread(dopplerSpreadHz, Fs, Nsam);
    end

    papr_log = [];
    for e=1:length(EbNodB)
        % generate a 2D array of QPSK symbols

        Nphases = 2^bps;
        tx_phases = pi/2*floor((rand(Nsym,Nc)*Nphases));
        if strcmp(method,"diversity")
          % duplicate carriers but with opposite phase
          tx_phases = [tx_phases (tx_phases-pi/2)];
        end
        tx_sym = gain*exp(j*tx_phases);

        % carrier frequencies, centre about 0
        st = floor(Nc*Nd/2);
        w = 2*pi/M*(-st:-st+Nc*Nd-1);
        
        tx = [];

        % generate OFDM signal

        for s=1:Nsym
          atx = zeros(1,M);
          for c=1:Nc*Nd
            atx += exp(j*(0:M-1)*w(c))*tx_sym(s,c);
          end
          % insert cyclic prefix and build up stream of time domain symbols
          tx = [tx atx(end-Ncp+1:end) atx];
        end
        Nsam = length(tx);

        if strcmp(channel,'multipath')
          assert(length(spread1) >= Nsam);
          assert(length(spread2) >= Nsam);
        end
        
        % bunch of PAPR reduction options
        
        tx_ = tx;
        if strcmp(method, "clip") || strcmp(method, "diversity")
          ind = find(abs(tx) > threshold);
          tx_(ind) = threshold*exp(j*angle(tx(ind)));
        end
        if strcmp(method, "compand1")
          tx_mag = interp1([0 1 3 5 6 20], [0 1 4 5.5 6 7], abs(tx), "pchip");
          tx_ = tx_mag.*exp(j*angle(tx));
        end
        if strcmp(method, "compand2")
          # power law compander x = a*y^power, y = (x/a) ^ (1/power)
          power=2; a=threshold/(threshold^power);
          tx_mag = (abs(tx)/a) .^ (1/power);
          tx_ = tx_mag.*exp(j*angle(tx));
        end
        
        if filt_en
          Nfilt=40;
          b = fir1(Nfilt,1.5*Nc*Nd/M);
          tx_ = filter(b,1,[tx_ zeros(1,Nfilt/2)]);
          tx_ = [tx_(Nfilt/2+1:end)];
        end

        rx = tx_;
        
        % multipath channel

        if strcmp(channel,'multipath')
          rx = spread1(1:Nsam).*rx + spread2(1:Nsam).*[zeros(1,path_delay) rx(path_delay+1:end)];
        end
        
        % normalise power after any multipath and non-linear shennanigans, so that noise addition is correct
        norm = sqrt(mean(abs(tx).^2)/mean(abs(rx).^2))
        rx *= norm;
        
        if phase_est
            % auxillary rx to get ideal phase ests on signal after multipath but before AWGN noise is added

            rx_phase = zeros(Nsym,Nc);
            for s=1:Nsym
              st = (s-1)*(M+Ncp)+1+timing; en = st+M-1;
              for c=1:Nc*Nd
                arx_sym = sum(exp(-j*(0:M-1)*w(c)) .* rx(st:en))/M;
                rx_phase(s,c) = arx_sym * conj(tx_sym(s,c));
              end
            end
            rx_phase = exp(j*arg(rx_phase));
        end

        % AWGN channel

        EsNodB = EbNodB(e) + 10*log10(bps);
        variance = M/(10^(EsNodB/10));
        noise = sqrt(variance/2)*randn(1,Nsam) + j*sqrt(variance/2)*randn(1,Nsam);
        rx += noise;

        % demodulate

        rx_sym = zeros(Nsym,Nc);
        for s=1:Nsym
          st = (s-1)*(M+Ncp)+1+timing; en = st+M-1;
          for c=1:Nc*Nd
            rx_sym(s,c) = sum(exp(-j*(0:M-1)*w(c)) .* rx(st:en))/M;
            if phase_est rx_sym(s,c) *= conj(rx_phase(s,c)); end
          end
          
          if strcmp(method,"diversity")
            for c=1:Nc
              rx_sym(s,c) += rx_sym(s,c+Nc)*exp(j*pi/2);
            end
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
          figure(2); clf; [hh nn] = hist(abs(tx),25,1);
          cdf = empirical_cdf((1:Nc),abs(tx));
          plotyy(nn,hh,1:Nc,cdf); title('PDF and CDF');
          % heat map type scatter plot
          figure(3); clf;
          rx_sym = reshape(rx_sym(:,1:Nc), Nsym*Nc, 1);
          h_data = [real(rx_sym)'; imag(rx_sym)']';
          h = log10(hist3(h_data,[50 50])+1);
          colormap("default"); imagesc(h);
          figure(4); clf; plot(real(rx_sym), imag(rx_sym), '+'); axis([-2 2 -2 2]);
          figure(5); clf; Tx_ = 10*log10(abs(fft(tx_))); plot(fftshift(Tx_));
          mx = 10*ceil(max(Tx_)/10); axis([1 length(Tx_) mx-60 mx]);
        end

        papr1 = 20*log10(max(abs(tx))/mean(abs(tx)));
        papr2 = 20*log10(max(abs(tx_))/mean(abs(tx_)));
        papr_log = [papr_log papr2];
        ber(e) = Terrs/Tbits;
        printf("EbNodB: %3.1f %3.1f PAPR: %5.2f %5.2f Tbits: %6d Terrs: %6d BER: %5.3f\n", EbNodB(e), 10*log10(mean(abs(tx_).^2)), papr1, papr2, Tbits, Terrs, ber(e))
    end

    papr = mean(papr_log);
end

% AWGN BER versus Eb/No curves -------------------------------------

function curves_awgn
    Nsym=1000;
    EbNodB=2:8;
    [ber1 papr1] = run_sim(Nsym,EbNodB, 0, filt_en=1);
    [ber2 papr2] = run_sim(Nsym,EbNodB, 0, filt_en=1, "clip", threshold=5);
    [ber3 papr3] = run_sim(Nsym,EbNodB, 0, filt_en=1, "clip", threshold=4);
    [ber4 papr4] = run_sim(Nsym,EbNodB, 0, filt_en=1, "compand1");
    [ber5 papr5] = run_sim(Nsym,EbNodB, 0, filt_en=1, "compand2", threshold=4);
    [ber6 papr6] = run_sim(Nsym,EbNodB, 0, filt_en=1, "diversity", threshold=3);

    figure(6); clf;
    semilogy(EbNodB, ber1,sprintf('b+-;vanilla OFDM %3.1f;',papr1),'markersize', 10, 'linewidth', 2); hold on;
    semilogy(EbNodB, ber2,sprintf('r+-;clip6 %3.1f;',papr2),'markersize', 10, 'linewidth', 2); 
    semilogy(EbNodB, ber3,sprintf('g+-;clip5 %3.1f;',papr3),'markersize', 10, 'linewidth', 2);
    semilogy(EbNodB, ber4,sprintf('c+-;compand1 %3.1f;',papr4),'markersize', 10, 'linewidth', 2);
    semilogy(EbNodB, ber5,sprintf('m+-;compand2 %3.1f;',papr5),'markersize', 10, 'linewidth', 2);
    semilogy(EbNodB, ber6,sprintf('bk+-;diversity %3.1f;',papr6),'markersize', 10, 'linewidth', 2);
    hold off;
    axis([min(EbNodB) max(EbNodB) 1E-4 1E-1]); grid;
    xlabel('Eb/No');
    print("papr_BER_EbNo.png","-dpng");

    figure(7); clf;
    semilogy(EbNodB+papr1, ber1,sprintf('b+-;vanilla OFDM %3.1f;',papr1),'markersize', 10, 'linewidth', 2); hold on;
    semilogy(EbNodB+papr2, ber2,sprintf('r+-;clip6 %3.1f;',papr2),'markersize', 10, 'linewidth', 2); 
    semilogy(EbNodB+papr3, ber3,sprintf('g+-;clip5 %3.1f;',papr3),'markersize', 10, 'linewidth', 2);
    semilogy(EbNodB+papr4, ber4,sprintf('c+-;compand1 %3.1f;',papr4),'markersize', 10, 'linewidth', 2);
    semilogy(EbNodB+papr5, ber5,sprintf('m+-;compand2 %3.1f;',papr5),'markersize', 10, 'linewidth', 2);
    semilogy(EbNodB+papr6, ber6,sprintf('bk+-;diversity %3.1f;',papr6),'markersize', 10, 'linewidth', 2);
    hold off;
    xlabel('Peak Eb/No');
    axis([6 20 1E-4 1E-1]); grid;
    print("papr_BER_peakEbNo.png","-dpng");
end


pkg load statistics;
more off;

% single point with lots of plots -----------

%run_sim(1000, EbNo=4, plot_en=1, filt_en=1, "diversity", threshold=3);
run_sim(1000, EbNo=10, channel='multipath', plot_en=1);

