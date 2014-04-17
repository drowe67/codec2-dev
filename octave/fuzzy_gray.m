% fuzzy_gray.m
% David Rowe
% 10 April 2014
%
% Experiments in fuzzy gray codes.  Idea is that with one bit error
% in the codeword only changes the encoded value by at most 1.

1;

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

function sim_out = test_baseline_uncoded(Ebvec, Nbits, Ntrials)
    Nlevels = 2.^ Nbits; powersOfTwo = 2 .^ fliplr(0:(Nbits-1));
    Nsymb   = Nbits;

    sim_out.qnoise_log = zeros(length(Ebvec),Ntrials);

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

function value = fuzzy_code_decode(codewords, rx_symbols)
    max_corr = 0;
    value = 0;
    [rows,cols] = size(codewords);
    for i=1:rows    
        corr = (-1 + 2*codewords(i,:)) * transpose(rx_symbols);
        if (corr > max_corr)
            max_corr = corr;
            value = i-1;
        end 
    end
endfunction

function valid_codewords = fuzzy_code_create(ndata, nparity)
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


function sim_out = test_fuzzy_code(Ebvec, Ndata, Ntrials)
    Nparity = 1;
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
        
        Terrs = 0;  Tbits = 0;
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

            rx_index = fuzzy_code_decode(codewords, rx_symbols);     
            rx_value = unquantise_value(rx_index, 0, 1, Nlevels);

            qsignal += tx_value*tx_value;
            qnoise  += (tx_value - rx_value) .^ 2;
            sim_out.qnoise_log(ne,nn) = tx_value - rx_value;
        end

        sim_out.QSNRvec(ne) = 10*log10(qsignal/qnoise);
        printf("EbNo (dB): %3.2f  QSNR (dB): %3.2f\n", EbNodB, sim_out.QSNRvec(ne));
    end
endfunction

Ebvec   = 0:7;
Ntrials = 1000;
Nbits   = 5;

baseline = test_baseline_uncoded(Ebvec, Nbits, Ntrials);
fuzzy = test_fuzzy_code(Ebvec, Nbits, Ntrials);

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
plot(baseline.qnoise_log(4,1:250),'b;baseline;')
hold on;
plot(fuzzy.qnoise_log(4,1:250),'r;fuzzy;')
hold off;
