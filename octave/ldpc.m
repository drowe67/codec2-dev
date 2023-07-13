% ldpc.m
%
#{
  David Rowe 2013
  Octave functions for the CML LDPC library.

  To install and compile CML support:
  
  $ sudo apt-get install liboctave-dev
  $ git clone git@github.com:drowe67/cml.git
  $ cd cml
  $ make

  If you have configured codec2 with cmake -DUNITTEST=1, then you will
  already have CML (e.g. under build_linux/cml), as it is used to run unit tests.

  To use CML when running Octave simulations from the Octave CLI, set an
  environment variable for CML_PATH in your shell or in your
  codec2/octave/.octaverc file:

    setenv("CML_PATH",sprintf("%s/codec2/build_linux/cml",getenv("HOME")))
#}

1;

function init_cml()
  currentdir = pwd;
  
  path_to_cml = getenv("CML_PATH");

  if exist(path_to_cml, 'dir') == 7
    cd(path_to_cml)
    CmlStartup      
    cd(currentdir); 
  else
    printf("\n---------------------------------------------------\n");
    printf("Can't start CML in path: %s\n", path_to_cml);
    printf("See CML_PATH instructions at top of this script (ldpc.m)\n");
    printf("-----------------------------------------------------\n\n");
    assert(0);
  end
end

% init using built in WiMax or DVSB2 code

function code_param = ldpc_init_builtin(code, rate, framesize, modulation, mod_order, mapping, constellation)
    if strcmp(code,'wimax')
      [code_param.H_rows, code_param.H_cols, code_param.P_matrix] = InitializeWiMaxLDPC( rate, framesize,  0 );
    end
    if strcmp(code,'dvbs2')
      [code_param.H_rows, code_param.H_cols, code_param.P_matrix] = InitializeDVBS2( rate, framesize);
    end
    if nargin == 7
      code_param.S_matrix = constellation;
    else
      if length(mapping) == 0
        code_param.S_matrix = CreateConstellation( modulation, mod_order);
      else
        code_param.S_matrix = CreateConstellation( modulation, mod_order, mapping );
      end
    end  
    code_param.bits_per_symbol = log2(mod_order);

    code_param.ldpc_data_bits_per_frame = length(code_param.H_cols) - length(code_param.P_matrix);
    code_param.ldpc_parity_bits_per_frame = framesize - code_param.ldpc_data_bits_per_frame;
    code_param.ldpc_coded_bits_per_frame = framesize;

    code_param.data_bits_per_frame  = code_param.ldpc_data_bits_per_frame;
    code_param.coded_bits_per_frame = code_param.ldpc_coded_bits_per_frame;
    code_param.coded_syms_per_frame = code_param.coded_bits_per_frame/code_param.bits_per_symbol;
endfunction


% init using user supplied code

function [code_param framesize rate] = ldpc_init_user(HRA, modulation, mod_order, mapping, constellation)
    [Nr Nc] = size(HRA);  
    rate = (Nc-Nr)/Nc;
    framesize = Nc;
    [H_rows, H_cols] = Mat2Hrows(HRA); 
    code_param.H_rows = H_rows; 
    code_param.H_cols = H_cols;
    code_param.P_matrix = [];
    if nargin == 5
      code_param.S_matrix = constellation;
    else
      if length(mapping) == 0
        code_param.S_matrix = CreateConstellation( modulation, mod_order);
      else
        code_param.S_matrix = CreateConstellation( modulation, mod_order, mapping );
      end
    end  
    code_param.bits_per_symbol = log2(mod_order);

    code_param.ldpc_data_bits_per_frame = length(code_param.H_cols) - length(code_param.P_matrix);
    code_param.ldpc_parity_bits_per_frame = framesize - code_param.ldpc_data_bits_per_frame;
    code_param.ldpc_coded_bits_per_frame = framesize;

    % these variables support 1's stuffing (not using all data bits to lower code rate)
    code_param.data_bits_per_frame  = code_param.ldpc_data_bits_per_frame;
    code_param.coded_bits_per_frame = code_param.ldpc_coded_bits_per_frame;
    code_param.coded_syms_per_frame = code_param.coded_bits_per_frame/code_param.bits_per_symbol;
endfunction


function [codeword s] = ldpc_enc(data, code_param)
  if code_param.data_bits_per_frame != code_param.ldpc_data_bits_per_frame
    % optionally lower the code rate by "1's stuffing" - setting Nunused data bits to 1
    Nunused = code_param.ldpc_data_bits_per_frame - code_param.data_bits_per_frame;
    codeword = LdpcEncode([data ones(1,Nunused)], code_param.H_rows, code_param.P_matrix);
    % remove unused data bits from codeword, as they are known to the receiver and don't need to be transmitted
    codeword = [ codeword(1:code_param.data_bits_per_frame) codeword(code_param.ldpc_data_bits_per_frame+1:end) ];
  else
    codeword = LdpcEncode( data, code_param.H_rows, code_param.P_matrix );
  end
  s = Modulate( codeword, code_param.S_matrix );
endfunction


function [detected_data paritychecks] = ldpc_dec(code_param, max_iterations, ...
                                                 demod_type, decoder_type, r, ...
                                                 EsNo, fading)
    % handle "1's stuffing" case where we don't use all data bits
    Nunused = code_param.ldpc_data_bits_per_frame - code_param.data_bits_per_frame;

    symbol_likelihood = Demod2D( r, code_param.S_matrix, EsNo, fading);
    
    % initialize the extrinsic decoder input

    input_somap_c = zeros(1, code_param.ldpc_coded_bits_per_frame - Nunused);
    bit_likelihood = Somap( symbol_likelihood, demod_type, input_somap_c );
    
    input_decoder_c = bit_likelihood(1:(code_param.ldpc_coded_bits_per_frame-Nunused));
    
    % insert "very likely" LLRs for unused data bits (in 1's stuffing case)
    input_decoder_c = [input_decoder_c(1:code_param.data_bits_per_frame) ...
                       100*ones(1,Nunused) ...
                       input_decoder_c(code_param.data_bits_per_frame+1:end)];
    
    [x_hat paritychecks] = MpDecode( -input_decoder_c, code_param.H_rows, code_param.H_cols, ...
                                      max_iterations, decoder_type, 1, 1);
    [mx mx_ind] = max(paritychecks);
    detected_data = x_hat(mx_ind,:);    
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
