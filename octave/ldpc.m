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

function dummy
endfunction

function sim_out = ldpc_proc(sim_in, resfile)

Eprob = sim_in.Eprob;  

framesize = sim_in.framesize;
rate      = sim_in.rate;
modulation = sim_in.modulation;
mod_order = sim_in.mod_order;
mapping   = sim_in.mapping;


Lim_Ferrs = sim_in.Lim_Ferrs;
Ntrials   = sim_in.Ntrials;
Esvec     = sim_in.Esvec;
deb       = sim_in.deb;


demod_type = 0;
decoder_type = 0;
max_iterations = 100;
code_param.bits_per_symbol = log2(mod_order);
bps = code_param.bits_per_symbol;

[code_param.H_rows, code_param.H_cols, code_param.P_matrix] = InitializeWiMaxLDPC( rate, sim_in.framesize,  0 );

code_param.data_bits_per_frame = length(code_param.H_cols) - length( code_param.P_matrix );

code_param.S_matrix = CreateConstellation( modulation, mod_order, mapping );

errfilename = '/home/david/codec2-dev/octave/mod_test_2000_poor_4dB.err';
fin = fopen(errfilename, "rb");
err = fread(fin,Inf, "short");
length(err)
Ntrials = floor(length(err)/framesize)

for ne = 1:length(Esvec)
    Es = Esvec(ne);
    EsNo = 10^(Es/10);

    
    Terrs = 0;  Tbits =0;  Ferrs =0;
    for nn = 1: Ntrials
        
        data = round( rand( 1, code_param.data_bits_per_frame ) );
        
        codeword = LdpcEncode( data, code_param.H_rows, code_param.P_matrix );
        code_param.code_bits_per_frame = length( codeword );
        Nsymb = code_param.code_bits_per_frame/bps;
       
        % modulate
        s = Modulate( codeword, code_param.S_matrix );
        code_param.symbols_per_frame = length( s );
     
        s = Modulate(codeword, code_param.S_matrix );
        code_param.symbols_per_frame = length( s );
        
      
        variance = 1/(2*EsNo);
        noise = sqrt(variance)*( randn(1,code_param.symbols_per_frame) + ...
            j*randn(1,code_param.symbols_per_frame) );
        a=ones(1, code_param.symbols_per_frame);  
        r = a.*s + noise;

        Nr = length(r); 
%        erasures = rand(1,Nr)<Eprob; 
%        r(erasures) = 0;            
        
        st = (nn-1)*Nr+1;
        en = nn*Nr;
        evec = err(st:en);
        r(find(evec)) = 0;

        symbol_likelihood = Demod2D( r, code_param.S_matrix, EsNo);
         
        % initialize the extrinsic decoder input
        input_somap_c = zeros(1, code_param.code_bits_per_frame );
        bit_likelihood = Somap( symbol_likelihood, demod_type, input_somap_c );
        
        input_decoder_c = bit_likelihood(1:code_param.code_bits_per_frame);
        
        x_hat= MpDecode( -input_decoder_c, code_param.H_rows, code_param.H_cols, ...
            max_iterations, decoder_type, 1, 1);
        detected_data = x_hat(max_iterations,:);
        error_positions = xor( detected_data(1:code_param.data_bits_per_frame), data );
        Nerrs = sum( error_positions);
        
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
    
    if length(resfile)>0
        save(resfile,  'sim_in',  'sim_out',  'cparams');
        disp(['Saved results to ' resfile '  at Es =' num2str(Es) 'dB']);
    end
  end
  
endfunction


