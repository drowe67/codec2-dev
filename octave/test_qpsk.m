% test_qpsk.m
% David Rowe Feb 2014
%
% Single sample per symbol QPSK modem simulation, based on code by Bill Cowley
% Generates curves BER versus E/No curves for different modems.  Design to
% test perform initial tests on coherent demodulation for HF channels without
% building a full blown modem.  Lacks filtering, timing estimation, frame sync.
% 

1;

% main test function 

function sim_out = ber_test(sim_in, modulation)

    framesize        = sim_in.framesize;
    Ntrials          = sim_in.Ntrials;
    Esvec            = sim_in.Esvec;
    phase_offset     = sim_in.phase_offset;
    phase_est        = sim_in.phase_est;
    w_offset         = sim_in.w_offset;
    plot_scatter     = sim_in.plot_scatter;
    Rs               = sim_in.Rs;
    hf_sim           = sim_in.hf_sim;
    hf_spread_hz     = sim_in.hf_spread_hz;
    hf_delay_samples = sim_in.hf_delay_samples;

    bps              = 2;
    Nsymb            = framesize/bps;
    prev_sym_tx      = qpsk_mod([0 0]);
    prev_sym_rx      = qpsk_mod([0 0]);
    rx_symb_log      = [];

    Np               = 5;
    r_delay_line     = zeros(1,Np);
    s_delay_line     = zeros(1,Np);
    hf_delay_line    = zeros(1,Nsymb+hf_delay_samples);
    spread_main_phi  = 0;
    spread_delay_phi = 0;
    spread_main_phi_log = [];

    % convert "spreading" samples from 1kHz carrier at Fs to complex baseband at Rs

    Fs = 8000; Fc = 1000;
    fspread = fopen("../unittest/sine1k_2Hz_spread.raw","rb");
    spread1k = fread(fspread, "int16")/10000;

    % down convert to complex baseband
    spreadbb = spread1k.*exp(-j*(2*pi*Fc/Fs)*(1:length(spread1k))');

    % remove -2000 Hz image
    b = fir1(50, 5/Fs);
    spreadlpf = filter(b,1,spreadbb);

    % decimate to symbol rate
    spread = spreadlpf(1:floor(Fs/Rs):length(spreadlpf));

    sc = 1;

    % design root nyquist (root raised cosine) filter and init tx and rx filter states

    alpha = 0.5; T=1/Fs; Nfiltsym=7; M=Fs/Rs;
    if floor(Fs/Rs) != Fs/Rs
        printf("oversampling ratio must be an integer\n");
        exit;
    end
    hrn = gen_rn_coeffs(alpha, T, Rs, Nfiltsym, M);
    Nfilter = length(hrn);
    tx_filter_memory = zeros(1, Nfilter);
    rx_filter_memory = zeros(1, Nfilter);
    s_delay_line_filt = zeros(1,Nfiltsym);
 
    wc = 2*pi*1500/Fs;

    for ne = 1:length(Esvec)
        Es = Esvec(ne);
        EsNo = 10^(Es/10);
    
        variance = Fs/(2*Rs*EsNo);
        Terrs = 0;  Tbits = 0;  Ferrs = 0;
        printf("EsNo (dB): %f EsNo: %f variance: %f\n", Es, EsNo, variance);

        tx_filt_log = [];
        rx_filt_log = [];
        rx_baseband_log = [];
        tx_baseband_log = [];
        noise_log = [];

        tx_phase = rx_phase = 0;

        for nn = 1: Ntrials
                  
           tx_bits = round( rand( 1, framesize ) );

           % modulate

            s = zeros(1, Nsymb);
            for i=1:Nsymb
                tx_symb = qpsk_mod(tx_bits(2*(i-1)+1:2*i));
                %printf("shift: %f prev_sym: %f  ", tx_symb, prev_sym_tx);
                if strcmp(modulation,'dqpsk')
                    tx_symb *= prev_sym_tx;
                    %printf("tx_symb: %f\n", tx_symb);
                    prev_sym_tx = tx_symb;
                end 
                s(i) = tx_symb;
            end

            s_ch = s;

            % root nyquist filter symbols

            for k=1:Nsymb

               % tx filter symbols

               tx_filt = zeros(1,M);

               % tx filter each symbol, generate M filtered output samples for each symbol.
               % Efficient polyphase filter techniques used as tx_filter_memory is sparse

               tx_filter_memory(Nfilter) = s_ch(k);

               for i=1:M
                   tx_filt(i) = M*tx_filter_memory(M:M:Nfilter) * hrn(M-i+1:M:Nfilter)';
               end
               tx_filter_memory(1:Nfilter-M) = tx_filter_memory(M+1:Nfilter);
               tx_filter_memory(Nfilter-M+1:Nfilter) = zeros(1,M);
               tx_filt_log = [tx_filt_log tx_filt];
               
               % AWGN noise and phase/freq offset channel simulation

               noise = sqrt(variance)*( randn(1,M) + j*randn(1,M) );
               noise_log = [noise_log noise];
               rx_baseband = tx_filt.*exp(j*phase_offset) + noise;
               phase_offset += w_offset;
               
               % rx filter symbol

               rx_filter_memory(Nfilter-M+1:Nfilter) = rx_baseband;
               rx_filt = rx_filter_memory * hrn';
               rx_filter_memory(1:Nfilter-M) = rx_filter_memory(1+M:Nfilter);
               rx_filt_log = [rx_filt_log rx_filt];

               % delay in tx data to compensate for filtering

               s_delay_line_filt(1:Nfiltsym-1) = s_delay_line_filt(2:Nfiltsym);
               s_delay_line_filt(Nfiltsym) = s(k);
               tx_bits(2*(k-1)+1:2*k) = qpsk_demod(s_delay_line_filt(1));
               s(k) = s_delay_line_filt(1);   % input to phase est later

               s_ch(k) = rx_filt;               
            end

            %noise = sqrt(variance)*( randn(1,Nsymb) + j*randn(1,Nsymb) );
            %s_ch = s_ch.*exp(j*phase_offset) + noise;
            %phase_offset += w_offset;

            % Channel simulation

            if hf_sim
                s_ch = s_ch.*spread(sc:sc+Nsymb-1)';
                sc += Nsymb;
            end

            % coherent demod phase estimation and correction
            % need sliding window
            % use known (pilot) symbols
            % start with all symbols pilots, then gradually decimate, e.g. 1 in 5 pilots

            if phase_est
                for i=1:Nsymb

                    % delay line for phase est window

                    r_delay_line(1:Np-1) = r_delay_line(2:Np);
                    r_delay_line(Np) = s_ch(i);

                    % delay in tx data to compensate data for phase est window

                    s_delay_line(1:Np-1) = s_delay_line(2:Np);
                    s_delay_line(Np) = s(i);
                    tx_bits(2*(i-1)+1:2*i) = qpsk_demod(s_delay_line(floor(Np/2)+1));

                    % estimate phase from known symbols and correct

                    corr = s_delay_line * r_delay_line';
                    s_ch(i) = r_delay_line(floor(Np/2)+1).*exp(j*angle(corr));
               end    
                %printf("corr: %f angle: %f\n", corr, angle(corr));
            end

            % de-modulate

            rx_bits = zeros(1, framesize);
            for i=1:Nsymb
                rx_symb = s_ch(i);
                if strcmp(modulation,'dqpsk')
                    tmp = rx_symb;
                    rx_symb *= conj(prev_sym_rx/abs(prev_sym_rx));
                    prev_sym_rx = tmp;
                end
                rx_bits((2*(i-1)+1):(2*i)) = qpsk_demod(rx_symb);
                rx_symb_log = [rx_symb_log rx_symb];
            end

            % Measure BER

            % discard bits from first 2*Nfiltsym symbols as tx and rx filter memories not full

            if nn == 1
                tx_bits = tx_bits(2*bps*Nfiltsym+1:length(tx_bits));
                rx_bits = rx_bits(2*bps*Nfiltsym+1:length(rx_bits));
            end

            error_positions = xor( rx_bits, tx_bits );
            Nerrs = sum(error_positions);
            Terrs += Nerrs;
            Tbits = Tbits + length(tx_bits);

        end
    
        TERvec(ne) = Terrs;
        FERvec(ne) = Ferrs;
        BERvec(ne) = Terrs/Tbits;
        printf("  Terrs: %d BER %f BER theory %f C %f N %f Es %f No %f Es/No %f\n\n", Terrs,
               Terrs/Tbits, 0.5*erfc(sqrt(EsNo/2)), var(tx_filt_log), var(noise_log),
               var(tx_filt_log)/Rs, var(noise_log)/Fs, (var(tx_filt_log)/Rs)/(var(noise_log)/Fs));
    end
    
    Ebvec = Esvec - 10*log10(bps);
    sim_out.BERvec = BERvec;
    sim_out.BER_theoryvec = 0.5*erfc(sqrt(10.^(Ebvec/10)));
    sim_out.Ebvec = Ebvec;
    sim_out.FERvec = FERvec;
    sim_out.TERvec  = TERvec;

    if plot_scatter
        figure(2);
        clf;
        scat = rx_symb_log(2*Nfiltsym:length(rx_symb_log)) .* exp(j*pi/4);
        plot(real(scat), imag(scat),'+');
    end
endfunction

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

% Start simulation ---------------------------------------

sim_in.Esvec            = 1:10; 
sim_in.Ntrials          = 100;
sim_in.framesize        = 30;
sim_in.phase_offset     = 0;
sim_in.phase_est        = 0;
sim_in.w_offset         = 0;
sim_in.plot_scatter     = 0;
sim_in.Rs               = 100;
sim_in.hf_sim           = 0;
sim_in.hf_spread_hz     = 2;
sim_in.hf_delay_samples = 5;

sim_qpsk                = ber_test(sim_in, 'qpsk');

sim_in.phase_offset     = 0;
sim_in.phase_est        = 0;
sim_in.w_offset         = 0;  
%sim_qpsk_coh            = ber_test(sim_in, 'qpsk');

sim_in.phase_offset     = 0;
sim_in.phase_est        = 0;
sim_in.w_offset         = 0;  
sim_in.plot_scatter     = 1;
sim_in.Esvec            = 7;
sim_in.hf_sim           = 0;
sim_qpsk_scatter        = ber_test(sim_in, 'qpsk');

figure(1); 
clf;
semilogy(sim_qpsk.Ebvec, sim_qpsk.BERvec)
hold on;
semilogy(sim_qpsk.Ebvec, sim_qpsk.BER_theoryvec,'r;coherent;')
%semilogy(sim_qpsk_coh.Ebvec, sim_qpsk_coh.BERvec,'r;coherent;')
hold off;
xlabel('Eb/N0')
ylabel('BER')
grid

