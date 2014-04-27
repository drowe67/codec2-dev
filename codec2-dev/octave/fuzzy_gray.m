% fuzzy_gray.m
% David Rowe
% 10 April 2014
%
% Various experiments in fuzzy gray codes and quantising and
% transmitting scalars. 

1;

% fuzzy gray coding idea: use an extra parity bit, if we get a single
% bit error the value will be "close: to the original, so effect of
% error will be soft.  Unlike data we don't need 0 bit errors. I
% struggled to extend this to larger m.

function three_bit_code
    m=4;
    log2_m=2;
    value_to_codeword = ["000"; "001"; "101"; "111"];
    codeword_to_value = [0 1 1 2 1 2 2 3 3];

    printf("tx_value tx_codeword rx_codeword rx_value distance\n");
    for i=1:m
        tx_codeword = bin2dec(value_to_codeword(i,:));
        tx_codeword_bin = value_to_codeword(i,:);
        rx_value = codeword_to_value(tx_codeword+1);
        distance = abs((i-1) - rx_value);
        printf("%8d %11s %11s %8d %8d\n", i-1, tx_codeword_bin, tx_codeword_bin, ...
               rx_value, distance );
    end
    printf("\n");
    for i=1:m
        tx_codeword = bin2dec(value_to_codeword(i,:));
        tx_codeword_bin = value_to_codeword(i,:);
        for j=1:(log2_m+1)
            rx_codeword = bitxor(tx_codeword, bitset(0,j));
            rx_codeword_bin = dec2bin(rx_codeword, 3);
            rx_value = codeword_to_value(rx_codeword+1);
            distance = abs((i-1) - rx_value);
            printf("%8d %11s %11s %8d %8d\n", i-1, tx_codeword_bin, rx_codeword_bin, ...
                   rx_value, distance );
       end
    end
endfunction

% regular natural binary quantiser

function index = quantise_value(value, min_value, max_value, num_levels)
    norm = (value - min_value)/(max_value - min_value);
    index = floor(num_levels * norm + 0.5);
    if (index < 0 ) 
        index = 0;
    end
    if (index > (num_levels-1)) 
        index = num_levels-1;
    end
endfunction

function value = unquantise_value(index, min_value, max_value, num_levels)
    step  = (max_value - min_value)/num_levels;
    value = min_value + step*(index);
endfunction

% converting natural binary to gray

function gray = binary_to_gray(natural)
    gray = bitxor(bitshift(natural,-1),natural);
endfunction

function natural = gray_to_binary(gray)
    for i=1:length(gray)
        mask = bitshift(gray(i),-1);
        num = gray(i);
        while(mask)
            num = bitxor(num, mask);
            mask = bitshift(mask,-1);
        end
        natural(i) = num;
    end
endfunction

function sim_out = test_baseline_uncoded(Ebvec, Nbits, Ntrials, enable_error_log, enable_gray)
    Nlevels = 2.^ Nbits; powersOfTwo = 2 .^ fliplr(0:(Nbits-1));
    Nsymb   = Nbits;

    sim_out.qnoise_log = zeros(length(Ebvec),Ntrials);
    sim_out.error_log  = [];

    for ne = 1:length(Ebvec)
        EbNodB = Ebvec(ne);
        EbNo = 10^(EbNodB/10);
    
        variance = 1/EbNo;
        
        Terrs = 0;  Tbits = 0;
        qsignal = qnoise = 0;
 
        for nn = 1:Ntrials
                  
            tx_value = rand(1,1);
            tx_index = quantise_value(tx_value, 0, 1, Nlevels);
            if enable_gray
                tx_index = binary_to_gray(tx_index);
            end
            tx_bits = dec2bin(tx_index, Nbits) - '0';
            tx_symbols = -1 + 2*tx_bits; 

            % AWGN noise and phase/freq offset channel simulation
            % 0.5 factor ensures var(noise) == variance , i.e. splits power between Re & Im

            noise = sqrt(variance*0.5)*(randn(1,Nsymb) + j*randn(1,Nsymb));
            rx_symbols = tx_symbols + noise;

            rx_bits = rx_symbols > 0;

            error_positions = xor(rx_bits, tx_bits);
            Nerrs = sum(error_positions);
            Terrs += Nerrs;
            Tbits += length(tx_bits);
            if enable_error_log
                sim_out.error_log  = [sim_out.error_log error_positions];
            end

            rx_index = (powersOfTwo  * rx_bits');
            if enable_gray
                rx_index = gray_to_binary(rx_index);
            end
            rx_value = unquantise_value(rx_index, 0, 1, Nlevels);

            qsignal += tx_value*tx_value;
            qnoise  += (tx_value - rx_value) .^ 2;
            sim_out.qnoise_log(ne,nn) = tx_value - rx_value;
        end

        sim_out.TERvec(ne) = Terrs;
        sim_out.BERvec(ne) = Terrs/Tbits;
        sim_out.QSNRvec(ne) = 10*log10(qsignal/qnoise);
        printf("EbNo (dB): %3.2f  Terrs: %6d BER %1.4f QSNR (dB): %3.2f\n", 
        EbNodB, Terrs, Terrs/Tbits, 10*log10(qsignal/qnoise));
    end

endfunction

function sim_out = test_varpower(Ebvec, Nbits, Ntrials, amps, enable_error_log)
    Nlevels = 2.^ Nbits; powersOfTwo = 2 .^ fliplr(0:(Nbits-1));
    Nsymb   = Nbits;

    sim_out.qnoise_log = zeros(length(Ebvec), Ntrials);
    sim_out.error_log  = [];

    for ne = 1:length(Ebvec)
        EbNodB = Ebvec(ne);
        EbNo = 10^(EbNodB/10);
    
        variance = 1/EbNo;
        
        Terrs = 0;  Tbits = 0;
        qsignal = qnoise = 0;
 
        for nn = 1:Ntrials
                  
            tx_value = rand(1,1);
            tx_index = quantise_value(tx_value, 0, 1, Nlevels);
            tx_bits = dec2bin(tx_index, Nbits) - '0';
            tx_symbols = (-1 + 2*tx_bits) .* amps;

            % AWGN noise and phase/freq offset channel simulation
            % 0.5 factor ensures var(noise) == variance , i.e. splits power between Re & Im

            noise = sqrt(variance*0.5)*(randn(1,Nsymb) + j*randn(1,Nsymb));
            rx_symbols = tx_symbols + noise;

            rx_bits = rx_symbols > 0;

            error_positions = xor(rx_bits, tx_bits);
            if enable_error_log
                sim_out.error_log  = [sim_out.error_log error_positions];
            end
            Nerrs = sum(error_positions);
            Terrs += Nerrs;
            Tbits += length(tx_bits);

            rx_index = (powersOfTwo  * rx_bits');
            rx_value = unquantise_value(rx_index, 0, 1, Nlevels);

            qsignal += tx_value*tx_value;
            qnoise  += (tx_value - rx_value) .^ 2;
            sim_out.qnoise_log(ne,nn) = tx_value - rx_value;
        end

        sim_out.TERvec(ne) = Terrs;
        sim_out.BERvec(ne) = Terrs/Tbits;
        sim_out.QSNRvec(ne) = 10*log10(qsignal/qnoise);
        printf("EbNo (dB): %3.2f  Terrs: %6d BER %1.4f QSNR (dB): %3.2f\n", 
        EbNodB, Terrs, Terrs/Tbits, 10*log10(qsignal/qnoise));
    end

endfunction

% gray codes with specified number of data and parity bits.  Soft
% decision decoding.  Didn't really work out.

function valid_codewords = fuzzy_code_create(ndata,nparity)
    Nbits = ndata + nparity;
    Nvalid = 2 .^ ndata;
    codewords = binary_to_gray(0:(2 .^ Nbits)-1);
    valid_codewords = dec2bin(codewords(1:2:(2 .^ Nbits)), Nbits) - '0';

    % check all valid codewords have a hamming distance of at least 2^nparity    

    bad_distance = 0;
    for i=1:Nvalid
        for k=i+1:Nvalid
            distance = sum(bitxor(valid_codewords(i,:), valid_codewords(k,:)));
            if distance < 2
                bad_distance++;
            end
        end
    end
    if bad_distance != 0
       printf("Error: Nvalid: %d bad_distance: %d\n", Nvalid, bad_distance);
       return;
    end

endfunction

function tx_codeword = fuzzy_code_encode(codewords, value)
    tx_codeword = codewords(value+1,:);
endfunction

function [value, best_codeword] = fuzzy_code_decode(codewords, rx_symbols)
    max_corr = 0;
    value = 0;
    [rows,cols] = size(codewords);
    for i=1:rows    
        corr = (-1 + 2*codewords(i,:)) * transpose(rx_symbols);
        if (corr > max_corr)
            max_corr = corr;
            value = i-1;
            best_codeword = codewords(i,:);
        end 
    end
endfunction


function sim_out = test_fuzzy_code(Ebvec, Ndata, Nparity, Ntrials)
    Nbits   = Ndata + Nparity;
    Nlevels = 2 .^ Ndata;
    Nsymb   = Nbits;
    powersOfTwo = 2 .^ fliplr(0:(Nbits-1));

    codewords = fuzzy_code_create(Ndata,Nparity);
    sim_out.qnoise_log = zeros(length(Ebvec), Ntrials);

    for ne = 1:length(Ebvec)
        EbNodB = Ebvec(ne);
        EbNo = 10^(EbNodB/10);
    
        variance = 1/EbNo;
        
        Terrs = 0; Terrs_coded = 0; Tbits = 0;
        Nsingle = Nsingle_corrected = 0;
        qsignal = qnoise = 0;
 
        for nn = 1:Ntrials
                  
            tx_value = rand(1,1);
            tx_index = quantise_value(tx_value, 0, 1, Nlevels);
            tx_codeword = fuzzy_code_encode(codewords, tx_index);
            tx_symbols = -1 + 2*tx_codeword; 

            % AWGN noise and phase/freq offset channel simulation
            % 0.5 factor ensures var(noise) == variance , i.e. splits power between Re & Im

            noise = sqrt(variance*0.5)*(randn(1,Nsymb) + j*randn(1,Nsymb));
            rx_symbols = tx_symbols + noise;

            % uncoded BER

            rx_bits = rx_symbols > 0;
            error_positions = xor(rx_bits(1:Ndata), tx_codeword(1:Ndata));
            Nerrs = sum(error_positions);
            Terrs += Nerrs;
            Tbits += Ndata;
            
            % decode and determine QSNR
 
            [rx_index, rx_codeword] = fuzzy_code_decode(codewords, rx_symbols);     
            rx_value = unquantise_value(rx_index, 0, 1, Nlevels);
            qsignal += tx_value*tx_value;
            qnoise  += (tx_value - rx_value) .^ 2;
            sim_out.qnoise_log(ne,nn) = tx_value - rx_value;

            % coded BER

            error_positions = xor(rx_codeword(1:Ndata), tx_codeword(1:Ndata));
            Nerrs_coded = sum(error_positions);
            if Nerrs == 1
                Nsingle++;
                if Nerrs_coded == 0
                    Nsingle_corrected++;
                end
            end
            Terrs_coded += Nerrs_coded;
 
        end

        sim_out.BERvec(ne) = Terrs/Tbits;
        sim_out.BERvec_coded(ne) = Terrs_coded/Tbits;
        sim_out.Nsingle(ne) = Nsingle;
        sim_out.Nsingle_corrected(ne) = Nsingle_corrected;

        sim_out.QSNRvec(ne) = 10*log10(qsignal/qnoise);
        printf("EbNo (dB): %3.2f  Terrs: %6d BER %1.4f Terrs_coded: %6d BER_coded %1.4f QSNR (dB): %3.2f", 
        EbNodB, Terrs, Terrs/Tbits, Terrs_coded, Terrs_coded/Tbits, sim_out.QSNRvec(ne));
        printf(" Nsingle: %d Nsingle_corrected: %d corrected: %3.1f\n", Nsingle,  Nsingle_corrected, Nsingle_corrected*100/Nsingle);
    end
endfunction

function compare_baseline_fuzzy
    Ebvec   = 0:3;
    Ntrials = 5000;
    Nbits   = 4; Nparity = 1;

    baseline = test_baseline_uncoded(Ebvec, Nbits, Ntrials, 0, 0);
    fuzzy = test_fuzzy_code(Ebvec, Nbits, Nparity, Ntrials);

    figure(1);
    clf;
    semilogy(Ebvec, baseline.BERvec)
    xlabel('Eb/N0')
    ylabel('BER')
    grid("minor")

    figure(2);
    clf;
    plot(Ebvec, baseline.QSNRvec,'b;baseline;')
    hold on;
    plot(Ebvec, fuzzy.QSNRvec,'r;fuzzy;')
    hold off;
    xlabel('Eb/N0')
    ylabel('SNR')
    grid("minor")

    figure(3);
    subplot(211)
    hist(baseline.qnoise_log(4,:),50);
    subplot(212)
    hist(fuzzy.qnoise_log(4,:),50);

    figure(4)
    subplot(211)
    plot(baseline.qnoise_log(4,1:250),'b;baseline;')
    subplot(212)
    plot(fuzzy.qnoise_log(4,1:250),'r;fuzzy;')
endfunction

% compare baseline and variable power schemes and make plots

function compare_baseline_varpower_plot
    Ebvec   = -2:5;
    Ntrials = 5000;
    Nbits   = 5;

    baseline = test_baseline_uncoded(Ebvec, Nbits, Ntrials, 0, 0);
    amps = [2 1.5 1.0 0.5 0.5];
    av_pwr = (amps*amps')/length(amps);
    amps_norm = amps/sqrt(av_pwr);
    varpower = test_varpower(Ebvec, Nbits, Ntrials,  amps_norm, 0);

    figure(1);
    clf;
    semilogy(Ebvec, baseline.BERvec)
    xlabel('Eb/No (dB)')
    ylabel('BER')
    grid("minor")
    title('BER versus Eb/No')

    figure(2);
    clf;
    plot(Ebvec, baseline.QSNRvec,'b;baseline;')
    hold on;
    plot(Ebvec, varpower.QSNRvec,'r;varpower;')
    hold off;
    xlabel('Eb/No (dB)')
    ylabel('SNR (dB)')
    grid("minor")
    title('Quantiser SNR versus Eb/No')

    figure(3);
    subplot(211)
    hist(baseline.qnoise_log(1,:),50);
    title('Baseline and Variable Power Error Histograms')
    subplot(212)
    hist(varpower.qnoise_log(1,:),50);

    figure(4)
    subplot(211)
    plot(baseline.qnoise_log(1,1:250),'b;baseline;')
    title('Baseline and Variable Power Error plots for Eb/No = -2dB')
    subplot(212)
    plot(varpower.qnoise_log(1,1:250),'r;varpower;')
endfunction

% Compare baseline and variable power schemes and make error pattern
% files for inserting into codec bit stream so we can listen to
% result.

function compare_baseline_varpower_error_files
    Ebvec   = -2;
    Fs      = 25;         % number of samples per second
    Nsec    = 15;         % seconds to simulate
    Ntrials = Fs*Nsec;
    Nbits   = 5;
    bits_per_frame = 52;
    bits_per_frame_rounded = ceil(bits_per_frame/8)*8; % c2enc uses integer number of bytes/frame
    start_bit = 12;                                    % first energy bit (after 4 voicing, 7 Wo bits)

    baseline = test_baseline_uncoded(Ebvec, Nbits, Ntrials, 1, 0);
    amps = [2 1.5 1.0 0.5 0.5];
    av_pwr = (amps*amps')/length(amps);
    amps_norm = amps/sqrt(av_pwr);
    varpower = test_varpower(Ebvec, Nbits, Ntrials, amps_norm, 1);

    % construct error patterns to apply to c2enc bit stream

    baseline_errors = [];
    for i=1:Ntrials
        error_positions = baseline.error_log(Nbits*(i-1)+1:Nbits*i);
        baseline_errors = [baseline_errors zeros(1,start_bit-1) error_positions ...
                           zeros(1, bits_per_frame_rounded - Nbits - (start_bit-1))];
    end
 
    varpower_errors = [];
    for i=1:Ntrials
        error_positions = varpower.error_log(Nbits*(i-1)+1:Nbits*i);
        varpower_errors = [varpower_errors zeros(1,start_bit-1) error_positions ...
                           zeros(1, bits_per_frame_rounded - Nbits - (start_bit-1))];
    end
 
    % save error patterns

    fep=fopen("energy_errors_baseline.bin","wb"); fwrite(fep, baseline_errors, "short"); fclose(fep);
    fep=fopen("energy_errors_varpower.bin","wb"); fwrite(fep, varpower_errors, "short"); fclose(fep);
endfunction

% compare natural and gray coding and make plots

function compare_natural_gray_plot
    Ebvec   = -2:10;
    Ntrials = 5000;
    Nbits   = 7;

    natural = test_baseline_uncoded(Ebvec, Nbits, Ntrials, 0, 0);
    gray    = test_baseline_uncoded(Ebvec, Nbits, Ntrials, 0, 1);

    figure(1);
    clf;
    semilogy(Ebvec, natural.BERvec)
    xlabel('Eb/No (dB)')
    ylabel('BER')
    grid("minor")
    title('BER versus Eb/No')

    figure(2);
    clf;
    plot(Ebvec, natural.QSNRvec,'b;natural;')
    hold on;
    plot(Ebvec, gray.QSNRvec,'r;gray;')
    hold off;
    xlabel('Eb/No (dB)')
    ylabel('SNR (dB)')
    grid("minor")
    title('Quantiser SNR versus Eb/No')

    figure(3);
    subplot(211)
    hist(natural.qnoise_log(1,:),50);
    title('Natural and Gray coded Error Histograms')
    subplot(212)
    hist(gray.qnoise_log(1,:),50);

    figure(4)
    subplot(211)
    plot(natural.qnoise_log(1,1:250),'b;natural;')
    axis([0 250 -1 1])
    title('Natural and Gray coded Error plots for Eb/No = -2dB')
    subplot(212)
    plot(gray.qnoise_log(1,1:250),'r;gray;')
    axis([0 250 -1 1])
endfunction

% compare natural at different Eb/No and Nbitsmake plots

function compare_natural_nbit_plot
    Ebvec   = -2:10;
    Ntrials = 5000;

    figure(1);
    clf;
    for n = 2:7
        natural = test_baseline_uncoded(Ebvec, n, Ntrials, 0, 0);
        plot(Ebvec, natural.QSNRvec)
        if n == 2
            hold on;
        end
   end
   hold off;

   xlabel('Eb/No (dB)')
   ylabel('SNR (dB)')
   grid("minor")
   title('Quantiser SNR versus Eb/No')
endfunction

function generate_varpower_error_files(EbNo, start_bit, end_bit, amps, error_file_name)
    Fs      = 25;         % number of samples per second
    Nsec    = 3;         % seconds to simulate
    Ntrials = Fs*Nsec;
    Nbits   = end_bit - start_bit + 1;
    bits_per_frame = 52;
    bits_per_frame_rounded = ceil(bits_per_frame/8)*8; % c2enc uses integer number of bytes/frame
                                                       % first energy bit (after 4 voicing, 7 Wo bits)

    % normalise powers and run test

    av_pwr = (amps*amps')/length(amps);
    amps_norm = amps/sqrt(av_pwr);
    av_pwr2 = (amps_norm*amps_norm')/length(amps_norm)
    varpower = test_varpower(EbNo, Nbits, Ntrials, amps_norm, 1);

    % construct error patterns to apply to c2enc bit stream

    varpower_errors = [];
    for i=1:Ntrials
        error_positions = varpower.error_log(Nbits*(i-1)+1:Nbits*i);

        if 0
        % reset single errors to tes effect of ideal single bit error correcting code
        for i=1:7
          st = 4*(i-1)+1
          en = 4*i
          if sum(error_positions(st:en)) == 1
            error_positions(st:en) = 0;
          end
        end
        for i=1:2
          st = 7*4+3*(i-1)+1 
          en = 7*4+3*i
          if sum(error_positions(st:en)) == 1
            error_positions(st:en) = 0;
          end
        end
        st = 7*4+3*2+1 
        en = 7*4+3*2+2
          if sum(error_positions(st:en)) == 1
            error_positions(st:en) = 0;
          end
        end

        num_errors(i) = sum(error_positions);
        varpower_errors = [varpower_errors zeros(1,start_bit-1) error_positions ...
                           zeros(1, bits_per_frame_rounded - Nbits - (start_bit-1))];
    end
 
    % save error pattern to file

    fep=fopen(error_file_name,"wb"); fwrite(fep, varpower_errors, "short"); fclose(fep);

    figure(1)
    clf
    hist(num_errors)
endfunction

more off;

%generate_varpower_error_files(0, 17, 52, ones(1,36), "lsp_baseline_errors_0dB.bin")
%amps = [1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 1 1 1 1 1 ];
%generate_varpower_error_files(0, 17, 52, amps, "lsp_varpower_errors_0dB.bin")

%compare_natural_nbit_plot
%compare_natural_gray_plot
%compare_baseline_varpower_plot
%compare_baseline_varpower_error_files

compare_baseline_fuzzy
%fuzzy_code_create(3,1)
