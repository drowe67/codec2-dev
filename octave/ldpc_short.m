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


function sim_out = run_sim(HRA, Ntrials, Esvec, genie_Es, packet_size, golay = 0)

  if golay
    rate = 0.5;
    code_param.data_bits_per_frame = 12;
    framesize = 24;
  else
    % LDPC code

    [Nr Nc] = size(HRA);  
    rate = (Nc-Nr)/Nc;
    framesize = Nc;
    modulation = 'BPSK';
    demod_type = 0;
    decoder_type = 0;
    max_iterations = 100;
    [H_rows, H_cols] = Mat2Hrows(HRA); 
    code_param.H_rows = H_rows; 
    code_param.H_cols = H_cols;
    code_param.P_matrix = [];
    code_param.data_bits_per_frame = length(code_param.H_cols) - length( code_param.P_matrix ); 
  end

  mod_order = 2; 
  bps = code_param.bits_per_symbol = log2(mod_order);

  % loop around each Esvec point
  
  for ne = 1:length(Esvec)
    Es = Esvec(ne);
    EsNo = 10^(Es/10);
    EbNodB = Es - 10*log10(code_param.bits_per_symbol * rate);

    Terrs = 0;  Tbits = 0;  Ferrs = 0; 
    tx_bits = rx_bits = [];

    for nn = 1: Ntrials        
      data = round( rand( 1, code_param.data_bits_per_frame ) );
      tx_bits = [tx_bits data];
      if golay
        codeword = egolayenc(data);
      else
        codeword = LdpcEncode( data, code_param.H_rows, code_param.P_matrix );
      end
      code_param.code_bits_per_frame = length( codeword );
      Nsymb = code_param.code_bits_per_frame/bps;      
       
      % modulate
      s = 1 - 2 * codeword;   
      code_param.symbols_per_frame = length( s );
              
      variance = 1/(2*EsNo);
      noise = sqrt(variance)* randn(1,code_param.symbols_per_frame); 
      r = s + noise;

      Nr = length(r);  

      if golay
        detected_data = egolaydec(r < 0);
        detected_data = detected_data(code_param.data_bits_per_frame+1:framesize);
      else

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
        detected_data = detected_data(1:code_param.data_bits_per_frame);
      end

      rx_bits = [rx_bits detected_data];

      error_positions = xor( detected_data, data );
      Nerrs = sum(error_positions);
        
      if Nerrs>0,  Ferrs = Ferrs +1;  end
      Terrs = Terrs + Nerrs;
      Tbits = Tbits + code_param.data_bits_per_frame;

    end
      
    % count packet errors using supplied packet size.  Operating one
    % one big string of bits lets us account for cases were FEC framesize
    % is less than packet size

    Perrs = 0; Tpackets = 0;
    for i=1:packet_size:length(tx_bits)-packet_size
      error_positions = xor( tx_bits(i:i+packet_size-1), rx_bits(i:i+packet_size-1));
      Nerrs = sum(error_positions);
      if Nerrs>0,  Perrs = Perrs +1;  end
      Tpackets++;
    end

    printf("EbNo: %3.1f dB BER: %5.4f PER: %5.4f Nbits: %4d Nerrs: %4d Tpackets: %4d Perr: %4d\n", EbNodB, Terrs/Tbits, Perrs/ Tpackets, Tbits, Terrs, Tpackets, Perrs);

    TERvec(ne) = Terrs;
    FERvec(ne) = Ferrs;
    BERvec(ne) = Terrs/ Tbits;
    PERvec(ne) = Perrs/ Tpackets;
    Ebvec = Esvec - 10*log10(code_param.bits_per_symbol * rate);
    
    sim_out.BERvec = BERvec;
    sim_out.PERvec = PERvec;
    sim_out.Ebvec = Ebvec;
    sim_out.FERvec = FERvec;
    sim_out.TERvec  = TERvec;
    sim_out.framesize = framesize;
  end
endfunction

% Start simulation here ----------------------------------------------

rand('seed',1);
randn('seed',1);
more off;
format;
init_cml;

Ntrials = 500;
Esvec = -3:0.5:3; 
packet_size = 28;


load HRA_112_112.txt
load HRA_112_56.txt
load HRA_56_56.txt
load HRA_56_28.txt

sim_out1 = run_sim(HRA_112_112, Ntrials, Esvec, 1, packet_size);
sim_out2 = run_sim(HRA_112_56 , Ntrials, Esvec, 1, packet_size);
sim_out3 = run_sim(HRA_56_56  , Ntrials*2, Esvec, 1, packet_size);
sim_out4 = run_sim(HRA_56_28  , Ntrials*2, Esvec, 1, packet_size);
sim_out5 = run_sim([], Ntrials*10, Esvec, 1, packet_size, 1);

Ebvec_theory = -2:0.5:6;
uncoded_BER_theory = 0.5*erfc(sqrt(10.^(Ebvec_theory/10)));

% need standard packet size to compare
% packet error if bit 0, or bit 1, or bit 2 .....
%              or bit 0 and bit 1
% no packet error if all bits ok (1-p(0))*(1-p(1))
% P(packet error) = p(0)+p(1)+....

uncoded_PER_theory = 1 - (1-uncoded_BER_theory).^packet_size;

figure(1); clf;
semilogy(Ebvec_theory,  uncoded_BER_theory, 'b+-;BPSK theory;','markersize', 10, 'linewidth', 2)
hold on;
semilogy(sim_out1.Ebvec, sim_out1.BERvec, 'g+-;rate 1/2 HRA 112 112;','markersize', 10, 'linewidth', 2)
semilogy(sim_out2.Ebvec, sim_out2.BERvec, 'r+-;rate 2/3 HRA 112 56;','markersize', 10, 'linewidth', 2)
semilogy(sim_out3.Ebvec, sim_out3.BERvec, 'c+-;rate 1/2 HRA 56 56;','markersize', 10, 'linewidth', 2)
semilogy(sim_out4.Ebvec, sim_out4.BERvec, 'k+-;rate 2/3 HRA 56 28;','markersize', 10, 'linewidth', 2)
semilogy(sim_out5.Ebvec, sim_out5.BERvec, 'm+-;rate 1/2 Golay (24,12);','markersize', 10, 'linewidth', 2)
hold off;
xlabel('Eb/No')
ylabel('BER')
grid
legend("boxoff");


figure(2); clf;
semilogy(Ebvec_theory,  uncoded_PER_theory, 'b+-;BPSK theory;','markersize', 10, 'linewidth', 2)
hold on;
semilogy(sim_out1.Ebvec, sim_out1.PERvec, 'g+-;rate 1/2 HRA 112 112;','markersize', 10, 'linewidth', 2)
semilogy(sim_out2.Ebvec, sim_out2.PERvec, 'r+-;rate 2/3 HRA 112 56;','markersize', 10, 'linewidth', 2)
semilogy(sim_out3.Ebvec, sim_out3.PERvec, 'c+-;rate 1/2 HRA 56 56;','markersize', 10, 'linewidth', 2)
semilogy(sim_out4.Ebvec, sim_out4.PERvec, 'k+-;rate 2/3 HRA 56 28;','markersize', 10, 'linewidth', 2)
semilogy(sim_out5.Ebvec, sim_out5.PERvec, 'm+-;rate 1/2 Golay (24,12);','markersize', 10, 'linewidth', 2)
hold off;
xlabel('Eb/No')
ylabel('PER')
grid
legend("boxoff");



