% lpdc_fsk_lib.m
% April 2015
%
% Library version of ldpc4.m written by vk5dsp.  Application is high bit rate
% balloon telemtry
%
% LDPC demo
% Call the CML routines and simulate one set of SNRs.  
% This fucntion is an updated version of ldpc3() which uses less 
% of the CML functions 
%
% sim_in the input parameter structure
% sim_out contains BERs and other stats for each value of SNR
% resfile is the result file
%

1;

function sim_out = ldpc5(sim_in, resfile, testmode, genie_Es, logging=0);

    if nargin<4, testmode = 0; end
    estEsN0 = 0; 

    HRA       = sim_in.HRA;
    framesize = sim_in.framesize;
    rate      = sim_in.rate;
    mod_order = sim_in.mod_order;

    Lim_Ferrs = sim_in.Lim_Ferrs;
    Ntrials   = sim_in.Ntrials;
    Esvec     = sim_in.Esvec;

    demod_type = 0;
    decoder_type = 0;
    max_iterations = 100;
    code_param = ldpc_init(HRA, mod_order);
    bps = code_param.bits_per_symbol;


    if (logging) 
       fod = fopen('decode.log', 'w'); 
       fwrite(fod, 'Es estEs Its  secs \n'); 
    end 


    for ne = 1:length(Esvec)
        Es = Esvec(ne);
        EsNo = 10^(Es/10);


        Terrs = 0;  Tbits =0;  Ferrs =0;
        for nn = 1: Ntrials

          data = round( rand( 1, code_param.data_bits_per_frame ) );
          codeword = ldpc_encode(code_param, data);

          code_param.code_bits_per_frame = length( codeword );
          Nsymb = code_param.code_bits_per_frame/bps;

            if testmode==1
               f1 = fopen("dat_in2064.txt", "w");
               for k=1:length(data);  fprintf(f1, "%u\n", data(k)); end
               fclose(f1); 
               system("./ra_enc"); 

               load("dat_op2064.txt"); 
               pbits = codeword(length(data)+1:end);   %   print these to compare with C code 
               dat_op2064(1:16)', pbits(1:16)  
               differences_in_parity =  sum(abs(pbits - dat_op2064'))
               pause;  
            end


            % modulate
            % s = Modulate( codeword, code_param.S_matrix );
            s= 1 - 2 * codeword;   
            code_param.symbols_per_frame = length( s );

            variance = 1/(2*EsNo);
            noise = sqrt(variance)* randn(1,code_param.symbols_per_frame); 
            %  +  j*randn(1,code_param.symbols_per_frame) );
            r = s + noise;
            Nr = length(r);  

            [detected_data Niters] = ldpc_decode(r, code_param, max_iterations, decoder_type);

            error_positions = xor( detected_data(1:code_param.data_bits_per_frame), data );
            Nerrs = sum( error_positions);

            t = clock;   t =  fix(t(5)*60+t(6)); 
            if (logging)  
              fprintf(fod, ' %3d  %4d\n', Niters, t); 
              end 

            if Nerrs>0, fprintf(1,'x'),  else fprintf(1,'.'),  end
            if (rem(nn, 50)==0),  fprintf(1,'\n'),  end

            if Nerrs>0,  Ferrs = Ferrs +1;  end
            Terrs = Terrs + Nerrs;
            Tbits = Tbits + code_param.data_bits_per_frame;

            if Ferrs > Lim_Ferrs, disp(['exit loop with #cw errors = ' ...
                    num2str(Ferrs)]);  break,  end
        end

        TERvec(ne) = Terrs;
        FERvec(ne) = Ferrs;
        BERvec(ne) = Terrs/ Tbits;
        Ebvec = Esvec - 10*log10(code_param.bits_per_symbol * rate);

        cparams= [code_param.data_bits_per_frame  code_param.symbols_per_frame ...
            code_param.code_bits_per_frame];

        sim_out.BERvec = BERvec;
        sim_out.Ebvec = Ebvec;
        sim_out.FERvec = FERvec;
        sim_out.TERvec  = TERvec;
        sim_out.cpumins = cputime/60;

        if nargin > 2
            save(resfile,  'sim_in',  'sim_out',  'cparams');
            disp(['Saved results to ' resfile '  at Es =' num2str(Es) 'dB']);
        end
      end
end


function code_param = ldpc_init(HRA, mod_order)
  code_param.bits_per_symbol = log2(mod_order);
  [H_rows, H_cols] = Mat2Hrows(HRA); 
  code_param.H_rows = H_rows; 
  code_param.H_cols = H_cols;
  code_param.P_matrix = [];
  code_param.data_bits_per_frame = length(code_param.H_cols) - length( code_param.P_matrix ); 
  code_param.symbols_per_frame = length(HRA);
end


function codeword = ldpc_encode(code_param, data)
      codeword = LdpcEncode( data, code_param.H_rows, code_param.P_matrix );
endfunction


% LDPC decoder - note it estimates EsNo from received symbols

function [detected_data Niters] = ldpc_decode(r, code_param, max_iterations, decoder_type)
  % in the binary case the LLRs are just a scaled version of the rx samples ..

  r = r / mean(abs(r));       % scale for signal unity signal  
  estvar = var(r-sign(r)); 
  estEsN0 = 1/(2* estvar + 1E-3); 
  input_decoder_c = 4 * estEsN0 * r;

  [x_hat, PCcnt] = MpDecode( input_decoder_c, code_param.H_rows, code_param.H_cols, ...
                             max_iterations, decoder_type, 1, 1);         
  Niters = sum(PCcnt!=0);
  detected_data = x_hat(Niters,:);

  if isfield(code_param, "c_include_file")
    write_code_to_C_include_file(code_param, max_iterations, decoder_type, input_decoder_c, x_hat);
  end
end


% optionally create a C include file for use in mpdecode.c C cmd line LDPC decoder

function write_code_to_C_include_file(code_param, max_iterations, decoder_type, input_decoder_c, x_hat)
       
  f = fopen(code_param.c_include_file, "wt");

  fprintf(f, "/*\n  FILE....: %s\n\n  Static arrays for TRC HF modem, generated", code_param.c_include_file);
  fprintf(f, "\n  by test_ldpc_fsk.m:simple_ut().\n\n*/\n\n");
   
  fprintf(f,"#define NUMBERPARITYBITS %d\n", rows(code_param.H_rows));
  fprintf(f,"#define MAX_ROW_WEIGHT %d\n", columns(code_param.H_rows));
  fprintf(f,"#define CODELENGTH %d\n", code_param.symbols_per_frame);
  fprintf(f,"#define NUMBERROWSHCOLS %d\n", rows(code_param.H_cols));
  fprintf(f,"#define MAX_COL_WEIGHT %d\n", columns(code_param.H_cols));
  fprintf(f,"#define DEC_TYPE %d\n", decoder_type);
  fprintf(f,"#define MAX_ITER %d\n", max_iterations);

  fprintf(f,"\ndouble H_rows[] = {\n");
     
  % clock out 2D array to linear C array in row order ....

  [r c] = size(code_param.H_rows);
  for j=1:c
    for i=1:r
      fprintf(f, "%f", code_param.H_rows(i,j));
      if (i == r) && (j ==c)
        fprintf(f,"\n};\n");
      else
        fprintf(f,", ");
      end
    end
  end

  fprintf(f,"\ndouble H_cols[] = {\n");
  [r c] = size(code_param.H_cols);
  for j=1:c
    for i=1:r
      fprintf(f, "%f", code_param.H_cols(i,j));
      if (i == r) && (j == c)
        fprintf(f,"\n};\n");
      else
        fprintf(f,", ");
      end
    end
  end

  fprintf(f,"\ndouble input[] = {\n");
  for i=1:length(input_decoder_c)
    fprintf(f, "%f", input_decoder_c(i));
    if i == length(input_decoder_c)
      fprintf(f,"\n};\n");
    else
      fprintf(f,", ");            
    end
  end

  fprintf(f,"\nint output[] = {\n");
  [r c] = size(x_hat);
  for j=1:c
    for i=1:r
      fprintf(f, "%d", x_hat(i,j));
      if (i == r) && (j == c)
        fprintf(f,"\n};\n");
      else
        fprintf(f,", ");
      end
    end
  end

  fclose(f);
end
