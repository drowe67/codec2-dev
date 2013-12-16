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

function sim_out = ldpc_proc(sim_in, resfile)

    rate = 3/4; 
    framesize = 576;  

    mod_order = 4; 
    modulation = 'QPSK';
    mapping = 'gray';

    demod_type = 0;
    decoder_type = 0;
    max_iterations = 100;

    code_param = ldpc_init(rate, framesize, modulation, mod_order, mapping);

    Ntrials = 84;
    EsNo=10;

    Tbits = Terrs = Ferrs = 0;
    
    data = [];
    r = []; 
    for nn = 1: Ntrials        
        d = round( rand( 1, code_param.data_bits_per_frame ) );
        data = [data d];
        [codeword, s] = ldpc_enc(d, code_param);
        code_param.code_bits_per_frame = length(codeword);
        code_param.symbols_per_frame = length(s);
        r = [r s];
    end

    for nn = 1: Ntrials        
        st = (nn-1)*code_param.symbols_per_frame + 1;
        en = (nn)*code_param.symbols_per_frame;
        detected_data = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, r(st:en), EsNo);
        st = (nn-1)*code_param.data_bits_per_frame + 1;
        en = (nn)*code_param.data_bits_per_frame;
        error_positions = xor( detected_data(1:code_param.data_bits_per_frame), data(st:en) );
        Nerrs = sum( error_positions);
        
        if Nerrs>0, fprintf(1,'x'),  else fprintf(1,'.'),  end
        if (rem(nn, 50)==0),  fprintf(1,'\n'),  end    
        if Nerrs>0,  Ferrs = Ferrs +1;  end
        Terrs = Terrs + Nerrs;
        Tbits = Tbits + code_param.data_bits_per_frame;        
    end

endfunction


