% ldpc.m
%
% David Rowe Mar 2017
% Based on Bill Cowley's LDPC simulations
%
% Octave simulation of BPSK with short LDPC codes developed by Bill.  First step
% in use of LDPC codes with FreeDV and Codec 2 700C.
%
% See lpdc.m for instruction son how to install CML library

1;

function init_cml
  currentdir = pwd;
  thiscomp = computer;

  if strfind(thiscomp, 'pc-linux-gnu')==8 
     if exist('LdpcEncode')==0, 
        %cd '~/cml'
        cd '/home/david/Desktop/cml'
        CmlStartup      
     end
  end
  cd(currentdir); 
end


function sim_out = ldpc4(HRA, Ntrials, Esvec, genie_Es);

  [Nr Nc] = size(HRA);  

  rate = (Nc-Nr)/Nc;
  framesize = Nc;
  mod_order = 2; 
  modulation = 'BPSK';
  demod_type = 0;
  decoder_type = 0;
  max_iterations = 100;
  bps = code_param.bits_per_symbol = log2(mod_order);

  [H_rows, H_cols] = Mat2Hrows(HRA); 
  code_param.H_rows = H_rows; 
  code_param.H_cols = H_cols;
  code_param.P_matrix = [];
  code_param.data_bits_per_frame = length(code_param.H_cols) - length( code_param.P_matrix ); 


  for ne = 1:length(Esvec)
    Es = Esvec(ne);
    EsNo = 10^(Es/10);
    EbNodB = Es - 10*log10(code_param.bits_per_symbol * rate);

    Terrs = 0;  Tbits = 0;  Ferrs = 0;

    for nn = 1: Ntrials        
      data = round( rand( 1, code_param.data_bits_per_frame ) );
      codeword = LdpcEncode( data, code_param.H_rows, code_param.P_matrix );
      code_param.code_bits_per_frame = length( codeword );
      Nsymb = code_param.code_bits_per_frame/bps;      
       
      % modulate
      s = 1 - 2 * codeword;   
      code_param.symbols_per_frame = length( s );
              
      variance = 1/(2*EsNo);
      noise = sqrt(variance)* randn(1,code_param.symbols_per_frame); 
      r = s + noise;

      Nr = length(r);  
      % in the binary case the LLRs are just a scaled version of the rx samples ...
      % if the EsNo is known -- use the following 
      if (genie_Es) 
	input_decoder_c = 4 * EsNo * r;   
      else
        r = r / mean(abs(r));       % scale for signal unity signal  
	estvar = var(r-sign(r)); 
	estEsN0 = 1/(2* estvar); 
	input_decoder_c = 4 * estEsN0 * r;
      end

      [x_hat, PCcnt] = MpDecode( input_decoder_c, code_param.H_rows, code_param.H_cols, ...
                                 max_iterations, decoder_type, 1, 1);
      Niters = sum(PCcnt!=0);
      detected_data = x_hat(Niters,:);
      error_positions = xor( detected_data(1:code_param.data_bits_per_frame), data );
      Nerrs = sum( error_positions);
        
      if Nerrs>0,  Ferrs = Ferrs +1;  end
      Terrs = Terrs + Nerrs;
      Tbits = Tbits + code_param.data_bits_per_frame;
    end
      
    printf("EbNo: %3.1f dB BER: %5.4f Nbits: %4d Nerrs: %4d\n", EbNodB, Terrs/Tbits, Tbits, Terrs);

    TERvec(ne) = Terrs;
    FERvec(ne) = Ferrs;
    BERvec(ne) = Terrs/ Tbits;
    Ebvec = Esvec - 10*log10(code_param.bits_per_symbol * rate);
    
    sim_out.BERvec = BERvec;
    sim_out.Ebvec = Ebvec;
    sim_out.FERvec = FERvec;
    sim_out.TERvec  = TERvec;
    sim_out.framesize = framesize;
  end
endfunction

% Start simulation here ----------------------------------------------

more off;
format;
init_cml

Ntrials =  2000;
Esvec = -3:0.5:3; 

sim_out1 = ldpc4(H1, Ntrials, Esvec, 1);
sim_out2 = ldpc4(H2, Ntrials, Esvec, 1);
sim_out3 = ldpc4(H3, Ntrials, Esvec, 1);
sim_out4 = ldpc4(H4, Ntrials, Esvec, 1);

EbNodBvec = sim_out1.Ebvec;
uncoded_BER_theory = 0.5*erfc(sqrt(10.^(EbNodBvec/10)));
uncoded_PER_theory = uncoded_BER_theory*sim_out1.framesize;

figure(1); clf;
semilogy(EbNodBvec,  uncoded_BER_theory, 'b')
hold on;
semilogy(EbNodBvec,  sim_out1.BERvec, 'g')
semilogy(EbNodBvec,  sim_out2.BERvec, 'r')
semilogy(EbNodBvec,  sim_out3.BERvec, 'c')
semilogy(EbNodBvec,  sim_out4.BERvec, 'k')
hold off;
xlabel('Eb/No')
ylabel('BER')
grid

figure(2); clf;
semilogy(EbNodBvec,  uncoded_PER_theory, 'b')
hold on;
semilogy(EbNodBvec,  sim_out1.FERvec/Ntrials, 'g')
semilogy(EbNodBvec,  sim_out2.FERvec/Ntrials, 'r')
semilogy(EbNodBvec,  sim_out3.FERvec/Ntrials, 'c')
semilogy(EbNodBvec,  sim_out4.FERvec/Ntrials, 'k')
hold off;
xlabel('Eb/No')
ylabel('PER')
grid
