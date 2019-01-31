% ldpc_gen_c_h_file.m
% David Rowe Sep 2015, B. Van Slyke 2019
%
% Create .c and h files for use in LDPC decoders
%
% NOTE:  You'll need to install the CML library as a number of functions involved
%        in LDPC use it.  See ldpc.m for instructions in installing the CML
%        library.
%
% Inputs:
% First parameter  - name of file with LDPC codes
% Second parameter - , defaults to 100
% Third parameter  - decoder_type, defaults to 0
%
% Output: Two files with the same filename as the LDPC input, but with .c and .h
% extensions.
 
function ldpc_gen_c_h_file(varargin)
    
    ldpc  % load ldpc functions
    ldpc_fsk_lib  % for ldpc_encode
    
    % Assuming cml has been installed in the users' home folder, which is the
    % default install location
    init_cml('~/cml/');
    
    if nargin == 0
        printf("Error - you must specify a file containing the LDPC codes (e.g. HRA_112_112.txt).\n");
        return;
    end
    loadStr = varargin{1};
  
 
    % The ldpc variable name may not be what we want for a file/variable names, but
    % the load filename will be, so use it.
    [~,ldpcArrayName,ext] = fileparts(loadStr);
    includeFileName = strcat(ldpcArrayName, '.h');
    sourceFileName = strcat(ldpcArrayName,  '.c');
   
    % Get the ext of the file first.  If it's a txt, then do what we
    % are doing.  If .mat, then just load, knowing the variable is HRA
    if strcmp(ext, '.mat') == 1
        load(loadStr);
    else
        % When calling 'load' this way, it returns a struct.  The code assumes the
        % struct has one element, and the one/first element is the array
        % to process
        tempStruct = load(loadStr);
        b = fieldnames(tempStruct);
        ldpcArrayName = b{1,1};
        % extract the array from the struct
        HRA = tempStruct.(ldpcArrayName);
    endif

    
    max_iterations = 100;
    decoder_type = 0;
 
    % user overloads 
    if nargin == 2
        max_iterations = varargin{2};
    end
 
    if nargin == 3
        decoder_type = varargin{3};
    end
 
   
    % the tests are performed using BPSK modulation, but in practice codes can be used
    % with other modulation, e.g. QPSK
    mod_order = 2; modulation = 'BPSK'; mapping = 'gray';
   
    [code_param framesize rate] = ldpc_init_user(HRA, modulation, mod_order, mapping);
  
    code_length = code_param.symbols_per_frame;
    code_length % reported as 112.
  
    % *********************  test for enc/dec
   [input_decoder_c, detected_data] = genInputOutputData(code_param, max_iterations, decoder_type);
     
   
    % First, create the H file
    f = fopen(includeFileName, "wt");
    printHeader(f, includeFileName, ldpcArrayName, mfilename());
 
    fprintf(f,"#define %s_NUMBERPARITYBITS %d\n", ldpcArrayName, rows(code_param.H_rows));
    fprintf(f,"#define %s_MAX_ROW_WEIGHT %d\n", ldpcArrayName, columns(code_param.H_rows));
    fprintf(f,"#define %s_CODELENGTH %d\n", ldpcArrayName, code_param.symbols_per_frame);
    fprintf(f,"#define %s_NUMBERROWSHCOLS %d\n", ldpcArrayName, rows(code_param.H_cols));
    fprintf(f,"#define %s_MAX_COL_WEIGHT %d\n", ldpcArrayName, columns(code_param.H_cols));
    fprintf(f,"#define %s_DEC_TYPE %d\n", ldpcArrayName, decoder_type);
    fprintf(f,"#define %s_MAX_ITER %d\n", ldpcArrayName, max_iterations);
    fprintf(f,"\n");
    fprintf(f,"extern const uint16_t %s_H_rows[];\n", ldpcArrayName);
    fprintf(f,"extern const uint16_t %s_H_cols[];\n", ldpcArrayName);
  
    fprintf(f,"extern const float %s_input[];\n", ldpcArrayName);
    fprintf(f,"extern const char %s_detected_data[];\n\n", ldpcArrayName);
  
    fclose(f);
  
    
    % Then, the C file
    f = fopen(sourceFileName, "wt");
    printHeader(f, sourceFileName, ldpcArrayName, mfilename());
    fprintf(f, "#include <stdint.h>\n");
    fprintf(f, "#include \"%s\"\n", includeFileName);  
    
    % clock out 2D array to linear C array in row order ....
    fprintf(f,"\nconst uint16_t %s_H_rows[] = {\n", ldpcArrayName);
    [r c] = size(code_param.H_rows);
    for j=1:c
        for i=1:r
            fprintf(f, "%d", code_param.H_rows(i,j));
            if (i == r) && (j ==c)  % weird, this does nothing
                fprintf(f,"\n};\n");
            else
                fprintf(f,", ");
            end
        end
    end
 
    fprintf(f,"\nconst uint16_t %s_H_cols[] = {\n", ldpcArrayName);
    [r c] = size(code_param.H_cols);
    for j=1:c
        for i=1:r
            fprintf(f, "%d", code_param.H_cols(i,j));
            if (i == r) && (j == c)
                fprintf(f,"\n};\n");
            else
                fprintf(f,", ");
            end
        end
    end
 
    % input and detected_data arrays
    fprintf(f,"const float %s_input[] = {\n", ldpcArrayName);
    for i=1:length(input_decoder_c)
        fprintf(f, "%.17g", input_decoder_c(i));
        if i == length(input_decoder_c)
            fprintf(f,"\n};\n");
        else
            fprintf(f,", ");         
        end
    end
    fprintf(f,"const char %s_detected_data[] = {\n", ldpcArrayName);
    for i=1:length(detected_data)
        fprintf(f, "%d", detected_data(i));
        if i == length(detected_data)
            fprintf(f,"\n};\n");
        else
            fprintf(f,", ");         
        end
    end
 
    fclose(f);
endfunction
 
function printHeader(f, includeFileName, ldpcArrayName, mFilename)
    fprintf(f, "/*\n  FILE....: %s\n\n", includeFileName);
    fprintf(f, "  Static arrays for LDPC codec %s, generated by %s.m.\n*/\n\n", ldpcArrayName, mFilename);
endfunction
 
function [input_decoder_c, detected_data] = genInputOutputData(code_param, max_iterations, decoder_type)
  
    % borrowed from test_ldpc_fsk_lib.m, simple_ut
    EsNodB = 3;
    data = round( rand( 1, code_param.data_bits_per_frame ) );
    codeword = ldpc_encode(code_param, data);    %defined in ldps_fsk_lib.
 
    s = 1 - 2 * codeword;  
    %aa = code_param.symbols_per_frame  % test - this WAS 112, then it gets set to 224???
    %code_param.symbols_per_frame = length( s )
 
    EsNo = 10^(EsNodB/10);
    variance = 1/(2*EsNo);
    noise = sqrt(variance)* randn(1,code_param.symbols_per_frame);
    r = s + noise;
   
    % borrowed from ldpc_fsk_lib.m, ldpc_decode
    llr = sd_to_llr(r);
 
    [x_hat, PCcnt] = MpDecode(llr, code_param.H_rows, code_param.H_cols, ...
                            max_iterations, decoder_type, 1, 1);        
    Niters = sum(PCcnt!=0);
    detected_data = x_hat(Niters,:);
 
    input_decoder_c = llr;
endfunction
