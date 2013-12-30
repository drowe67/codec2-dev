% ldpc.m
% ldpc functions

% LDPC demo;   Bill Cowley 
% Call the CML routines and simulate one set of SNRs.   See test_ldpc1.m
%
% sim_in the input parameter structure
% sim_out contains BERs and other stats for each value of SNR
% resfile is the result file
%
% 4/oct/2013:   edited to use the WiMax eIRA codes 
%               see 'help InitializeWiMaxLDPC' for parameter values 

1;

function code_param = ldpc_init(rate, framesize, modulation, mod_order, mapping)
    [code_param.H_rows, code_param.H_cols, code_param.P_matrix] = InitializeWiMaxLDPC( rate, framesize,  0 );
    code_param.data_bits_per_frame = length(code_param.H_cols) - length( code_param.P_matrix );
    code_param.S_matrix = CreateConstellation( modulation, mod_order, mapping );
    code_param.bits_per_symbol = log2(mod_order);
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

% inserts a unique word into a frame of bits.  The UW bits are spread
% throughout the input frame 2 bits at a time.

function frameout = insert_uw(framein, uw)

    luw = length(uw);
    lframein = length(framein);
    spacing = 2*lframein/luw;

    frameout = [];
    mod_uw = [];

    pin = 1; pout = 1; puw = 1;
    while (luw)
        %printf("pin %d pout %d puw %d luw %d\n", pin, pout, puw, luw);
        frameout(pout:pout+spacing-1) = framein(pin:pin+spacing-1);
        pin += spacing; 
        pout += spacing;
        frameout(pout:pout+1) = uw(puw:puw+1);
        mod_uw(pout:pout+1) = qpsk_mod(uw(puw:puw+1));
        puw += 2;
        pout += 2;
        luw -= 2;
    end
endfunction

% removes a unique word from a frame of bits.  The UW bits are spread
% throughout the input frame 2 bits at a time.

function frameout = remove_uw(framein, lvd, luw)

    spacing = 2*lvd/luw;

    frameout = [];

    pin = 1; pout = 1;
    while (luw)
        %printf("pin %d pout %d luw %d  ", pin, pout, luw);
        %printf("pin+spacing-1 %d lvd %d lframein: %d\n", pin+spacing-1, lvd, length(framein));
        frameout(pout:pout+spacing-1) = framein(pin:pin+spacing-1);
        pin  += spacing + 2; 
        pout += spacing;
        luw  -= 2;
    end

endfunction


% builds up a sparse QPSK modulated version version of the UW for use
% in UW sync at the rx

function mod_uw = build_mod_uw(uw, spacing)
    luw = length(uw);

    mod_uw = [];

    pout = 1; puw = 1;
    while (luw)
        %printf("pin %d pout %d puw %d luw %d\n", pin, pout, puw, luw);
        pout += spacing/2;
        mod_uw(pout) = qpsk_mod(uw(puw:puw+1));
        puw += 2;
        pout += 1;
        luw -= 2;
    end
endfunction


% Uses the UW to determine when we have a full codeword ready for decoding

function found_uw = look_for_uw(mem_rx_symbols, mod_uw)
    sparse_mem_rx_symbols = mem_rx_symbols(find(mod_uw));

    % correlate with ref UW

    num = (mem_rx_symbols * mod_uw') .^ 2;
    den = (sparse_mem_rx_symbols * sparse_mem_rx_symbols') * (mod_uw * mod_uw');

    found_uw = abs(num/(den+1E-6)) > 0.8;
endfunction


function [codeword s] = ldpc_enc(data, code_param)
        codeword = LdpcEncode( data, code_param.H_rows, code_param.P_matrix );
        s = Modulate( codeword, code_param.S_matrix );
endfunction


function detected_data = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, r, EsNo)
    symbol_likelihood = Demod2D( r, code_param.S_matrix, EsNo);
         
    % initialize the extrinsic decoder input
    input_somap_c = zeros(1, code_param.code_bits_per_frame );
    bit_likelihood = Somap( symbol_likelihood, demod_type, input_somap_c );
        
    input_decoder_c = bit_likelihood(1:code_param.code_bits_per_frame);
        
    x_hat= MpDecode( -input_decoder_c, code_param.H_rows, code_param.H_cols, ...
                     max_iterations, decoder_type, 1, 1);
    detected_data = x_hat(max_iterations,:);
endfunction


% Packs a binary array into an array of 8 bit bytes, MSB first

function packed = packmsb(unpacked)
    packed = zeros(1,floor(length(unpacked)+7)/8);
    bit = 7; byte = 1;
    for i=1:length(unpacked)
        packed(byte) = bitor(packed(byte), bitshift(unpacked(i),bit));
        bit--;
        if (bit < 0)
            bit = 7;
            byte++;
        end 
    end
endfunction


% unpacks an array of 8 bit bytes into a binary array of unpacked bits, MSB first

function unpacked = unpackmsb(packed)
    bit = 7; byte = 1;
    for i=1:length(packed)*8
        unpacked(i) = bitand(bitshift(packed(byte), -bit), 1);
        bit--;
        if (bit < 0)
            bit = 7;
            byte++;
        end 
    end
endfunction

