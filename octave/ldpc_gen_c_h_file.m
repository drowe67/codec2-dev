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
% Second parameter - max_iterations, defaults to 100
% Third parameter  - decoder_type, defaults to 0
% Fourth parameter - input array, optional
% Fifth parameter  - detected data array, optional
%
% Output: Two files with the same filename as the LDPC input, but with .c and .h
% extensions.

function ldpc_gen_c_h_file(varargin)
 
    ldpc  % load ldpc functions
 
    if nargin == 0
        printf("Error - you must specify a file containing the LDPC codes (e.g. HRA_112_112.txt).\n");
        return;
    end
 
    loadStr = varargin{1};
   
    if nargin == 1
        max_iterations = 100;
        decoder_type = 0;
    else
        max_iterations = varargin{2};
        decoder_type = varargin{3};
    end
  
    mod_order = 4; bps = 2; modulation = 'QPSK'; mapping = 'gray';
    % Assuming cml has been installed in the users' home folder, which is the
    % default install location
    init_cml('~/cml/');
    
    % When calling 'load' this way, it returns a struct.  The code assumes the
    % struct has one element, and the one/first element is the array
    % to process
    tempStruct = load(loadStr);
    b = fieldnames(tempStruct);
    ldpcArrayName = b{1,1};
 
    % extract the array from the struct
    ldpcArray = tempStruct.(ldpcArrayName);
   
    % The ldpc variable name may not be what we want for a file/variable names, but
    % the load filename will be, so use it.
    [~,ldpcArrayName,~] = fileparts(loadStr);
    includeFileName = strcat(ldpcArrayName, '.h');
    sourceFileName = strcat(ldpcArrayName,  '.c');
    
    [code_param framesize rate] = ldpc_init_user(ldpcArray, modulation, mod_order, mapping);
   
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
   
    if nargin == 5
        fprintf(f,"extern const float %s_input[];\n", ldpcArrayName);
        fprintf(f,"extern const char %s_detected_data[];\n\n", ldpcArrayName);
    end
   
    fclose(f);
   
    
    % Then, the C file
    f = fopen(sourceFileName, "wt");
 
    printHeader(f, sourceFileName, ldpcArrayName, mfilename());
 
    fprintf(f, "#include <stdint.h>\n");
    fprintf(f, "#include \"%s\"\n", includeFileName);   
    
    % clock out 2D array to linear C array in row order ....
    fprintf(f,"\ndouble %s_H_rows[] = {\n", ldpcArrayName);
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
    fprintf(f,"\ndouble %s_H_cols[] = {\n", ldpcArrayName);
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
    % If the input and detected output have been specified, then print them
    if nargin == 5
      input_decoder_c = varargin{4};
      detected_data = varargin{5};
 
        fprintf(f,"\ndouble input[] = {\n");
        for i=1:length(input_decoder_c)
            fprintf(f, "%.17g", input_decoder_c(i));
            if i == length(input_decoder_c)
                fprintf(f,"\n};\n");
            else
                fprintf(f,", ");          
            end
        end
 
        fprintf(f,"\nchar detected_data[] = {\n");
        for i=1:length(detected_data)
            fprintf(f, "%d", detected_data(i));
            if i == length(detected_data)
                fprintf(f,"\n};\n");
            else
                fprintf(f,", ");          
            end
        end
    end % if nargin == 5
    fclose(f);
endfunction
 
 
function printHeader(f, includeFileName, ldpcArrayName, mFilename)
    fprintf(f, "/*\n  FILE....: %s\n\n", includeFileName);
    fprintf(f, "  Static arrays for LDPC codec %s, generated by %s.m.\n*/\n\n", ldpcArrayName, mFilename);
endfunction
