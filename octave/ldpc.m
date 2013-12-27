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

% inserts a unique word into a frame of bits

function frameout = insert_uw(framein, uw)

    luw = length(uw);
    lframein = length(framein);
    spacing = lframein/luw;

    frameout = [];

    for i=1:luw
        frameout(1+(i-1)*spacing+i-1:i*spacing+i-1) = framein(1+(i-1)*spacing:i*spacing);
        frameout(i*spacing+i) = uw(i);
    end

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

