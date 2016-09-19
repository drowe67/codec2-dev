% ldpc.m
%
% David Rowe 2013
% Octave functions to help us use the CML LDPC code.
%
% Installing CML library
% ----------------------
%
% $ sudo apt-get install liboctave-dev
% $ wget http://www.iterativesolutions.com/user/image/cml.1.10.zip
% $ unzip cml.1.10.zip
% $ patch -p0 < ~/codec2-dev/octave/cml.patch
% $ cd cml/source
% $ octave
% octave:> make
% (you'll see a few warnings but hopefully no errors)

1;


function code_param = ldpc_init(rate, framesize, modulation, mod_order, mapping)
    [code_param.H_rows, code_param.H_cols, code_param.P_matrix] = InitializeWiMaxLDPC( rate, framesize,  0 );
    code_param.data_bits_per_frame = length(code_param.H_cols) - length( code_param.P_matrix );
    code_param.S_matrix = CreateConstellation( modulation, mod_order, mapping );
    code_param.bits_per_symbol = log2(mod_order);
    code_param.code_bits_per_frame = framesize;
    code_param.symbols_per_frame = framesize/2;
endfunction



function [codeword s] = ldpc_enc(data, code_param)
        codeword = LdpcEncode( data, code_param.H_rows, code_param.P_matrix );
        s = Modulate( codeword, code_param.S_matrix );
endfunction


function detected_data = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, r, EsNo, fading)
    if nargin == 6
      fading = ones(1, length(r));
    end

    symbol_likelihood = Demod2D( r, code_param.S_matrix, EsNo, fading);
         
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


% symbol interleaver that acts on bits 2 at a time

function y = interleave_bits(interleaver, x)
    y =  zeros(1,length(x));
    for i = 1:length(interleaver)
        dst = interleaver(i);
        y(2*(dst-1)+1:2*dst) = x(2*(i-1)+1:2*(i));
    end
endfunction

% symbol de-interleaver

function x = deinterleave_symbols(interleaver, y)
    for i = 1:length(interleaver)
        x(i) = y(interleaver(i));
    end
endfunction
