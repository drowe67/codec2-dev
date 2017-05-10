% ldpc_short.m
%
% David Rowe Mar 2017
% Based on Bill Cowley's LDPC simulations
%
% Octave simulation of BPSK with short LDPC codes developed by Bill.  First step
% in use of LDPC codes with FreeDV and Codec 2 700C.  BPSK makes life a bits easier
% for simulation
%
% NOTE: You will need to set the CML path in the call to init_cml() below
%       for you CML install.  See lpdc.m for instructions on how to install 
%       CML library

1;

gp_interleaver;

function init_cml(path_to_cml)
  currentdir = pwd;
  
  if exist(path_to_cml, 'dir') == 7
    cd(path_to_cml)
    CmlStartup      
    cd(currentdir); 
  else
    printf("\n---------------------------------------------------\n");
    printf("Can't start CML in path: %s\n", path_to_cml);
    printf("See CML path instructions at top of this script\n");
    printf("-----------------------------------------------------\n\n");
    assert(0);
  end
end


function sim_out = run_sim(sim_in, HRA, Ntrials)

  rand('seed',1);
  randn('seed',1);

  % Note this is effective Eb/No of payload data bits, sorta thing we
  % plot on BER versus Eb/No graphs of decoded data.  So if we have a
  % rate 1/2 code, each codeword bit will have Eb/No - 3dB.

  EbNovec     = sim_in.EbNovec;

  genie_Es    = sim_in.genie_Es;
  code        = sim_in.code;
  hf_en       = sim_in.hf_en;
  verbose     = sim_in.verbose;

  if isfield(sim_in, "interleave_en") 
    interleave_en = sim_in.interleave_en;
    interleave_frames = sim_in.interleave_frames;
  else
    interleave_en = 0;
    interleave_frames = 1;
  end

  % default number of carriers for HF model, we fiddle with this a bit
  % for different FEC schemes

  Nc = 28;

  % set up for different FEC, get roughly the same frame size so for
  % HF channels the simulation runs over roughly the same time

  if strcmp(code, 'ldpc')
    [numr numc] = size(HRA);
    rate = (numc-numr)/numc;
    modulation = 'BPSK';
    demod_type = 0;
    decoder_type = 0;
    max_iterations = 100;
    [H_rows, H_cols] = Mat2Hrows(HRA); 
    code_param.H_rows = H_rows; 
    code_param.H_cols = H_cols;
    code_param.P_matrix = [];
    code_param.data_bits_per_frame = length(code_param.H_cols) - length( code_param.P_matrix ); 
    code_param.code_bits_per_frame = numc;
    Ncodewordsperframe = interleave_frames*224/numc;
    framesize = Ncodewordsperframe*numc;
  end

  if strcmp(code, 'golay')
    rate = 0.5;
    code_param.data_bits_per_frame = 12;
    code_param.code_bits_per_frame = 24;

    % one Golay (24,12) codeword per row

    Nc = 24;   
    Ncodewordsperframe = 9;
    framesize = Nc*Ncodewordsperframe;
  end

  if strcmp(code, 'diversity')
    rate = 0.5;
    code_param.data_bits_per_frame = Nc/2;
    code_param.code_bits_per_frame = Nc;
    Ncodewordsperframe = 8;
    framesize = Nc*Ncodewordsperframe;
  end

  Rs = 50; Tp = (framesize/Nc)/Rs; Tp_codec2 = 0.04;

  % we are using BPSK here

  mod_order = 2; 
  bps = code_param.bits_per_symbol = log2(mod_order);

  % init HF model

  if hf_en

    % some typical values, assume symbols spread across Nc OFDM carriers

    dopplerSpreadHz = 1; path_delay = 1E-3*Rs;

    nsymb = Ntrials*framesize*bps;
    spread1 = doppler_spread(dopplerSpreadHz, Rs, 1.1*ceil(nsymb/Nc));
    spread2 = doppler_spread(dopplerSpreadHz, Rs, 1.1*ceil(nsymb/Nc));
    hf_gain = 1.0/sqrt(var(spread1)+var(spread2));
    % printf("nsymb: %d lspread1: %d\n", nsymb, length(spread1));
    hf_model = zeros(ceil(nsymb/Nc),Nc);
  end

  % --------------------------------------------------------------------
  % Sun a simulation for every EbNovec point
  %---------------------------------------------------------------------

  for ne = 1:length(EbNovec)

    % Given Eb/No of payload data bits, work out Es/No we need to
    % apply to each codeword symbol, e.g. for rate 1/2 code:
    %
    % Es/No = Eb/No - 3 dB

    EbNodB = EbNovec(ne);
    EsNodB = EbNodB + 10*log10(code_param.bits_per_symbol*rate);
    EsNo = 10^(EsNodB/10);

    % reset Hf channel index here so each run gets exactly the same HF channel

    hf_r = 1;

    % bunch of counters and logs

    Terrs = 0;  Tbits = 0;  Ferrs = 0;  Terrs_raw = Tbits_raw = 0;
    r_log = []; Nerrs_log = []; Nerrs_raw_log = [];
    tx_bits_log = rx_bits_log = [];

    for nn = 1: Ntrials        

      % Each trial is one complete frame - OK first set up frame to tx

      codeword = tx_bits = [];

      if strcmp(code, 'ldpc')
        for f=1:Ncodewordsperframe
          data = round( rand( 1, code_param.data_bits_per_frame ) );
          tx_bits = [tx_bits data];
          codeword = [codeword LdpcEncode(data, code_param.H_rows, code_param.P_matrix) ];
        end
        if interleave_en
          codeword = gp_interleave(codeword);
        end
      end

      if strcmp(code, 'golay')
        for f=1:Ncodewordsperframe
          data = round( rand( 1, code_param.data_bits_per_frame ) );
          tx_bits = [tx_bits data];
          codeword = [codeword egolayenc(data)];
        end
      end

      if strcmp(code, 'diversity')
        for f=1:Ncodewordsperframe
          data = round( rand( 1, code_param.data_bits_per_frame ) );
          tx_bits = [tx_bits data];
          codeword = [codeword data data];
        end
      end

      tx_bits_log = [tx_bits_log tx_bits];

      Nsymb = code_param.code_bits_per_frame/bps;      
       
      % modulate

      s = 1 - 2 * codeword;   
      Ns = code_param.symbols_per_frame = length(s);
              
      if hf_en

        % Simplified rate Rs OFDM simulation model that doesn't
        % include ISI, just freq filtering.  We assume perfect phase
        % estimation so it's just amplitude distortion.  We distribute
        % symbols across Nc carriers

        Nr = ceil(length(s)/Nc);
        w = (1:Nc)*2*pi;  
        s = [s zeros(1,Nr*Nc-Ns)];  % pad out to integer number of rows

        for r=1:Nr
          for c=1:Nc
            hf_model(hf_r,c) = hf_gain*(spread1(hf_r) + exp(-j*w(c)*path_delay)*spread2(hf_r));
            s(Nc*(r-1)+c) *= abs(hf_model(hf_r,c));
          end
          hf_r++;
        end
        s = s(1:Ns); % remove padding
      end
      
      variance = 1/(EsNo);
      noise = sqrt(variance/2)*(randn(1,Ns) + j*randn(1,Ns));
      r = s + noise;

      Nr = length(r);  
      r_log = [r_log r];

      % raw BER

      detected_data = real(r) < 0;
      error_positions = xor(detected_data, codeword);

      Nerrs_raw = sum(error_positions);
      Terrs_raw += Nerrs_raw;
      Tbits_raw += length(codeword);
      Nerrs_raw_log = [Nerrs_raw_log Nerrs_raw];

      % FEC decoder

      rx_bits = [];

      if strcmp(code, 'ldpc')
        if interleave_en
          r = gp_deinterleave(r);
        end

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

        for f=1:Ncodewordsperframe
          st = (f-1)*code_param.code_bits_per_frame + 1;
          en = st + code_param.code_bits_per_frame - 1;
          [x_hat, PCcnt] = MpDecode(input_decoder_c(st:en), code_param.H_rows, code_param.H_cols, ...
                                   max_iterations, decoder_type, 1, 1);
          Niters = sum(PCcnt!=0);
          arx_bits = x_hat(Niters,:);
          rx_bits = [rx_bits arx_bits(1:code_param.data_bits_per_frame)];

          #{
            if isfield(sim_in, "c_include_file") 

            % optionally dump code and unit test data to a C header file

            code_param.c_include_file = sim_in.c_include_file;
            ldpc_gen_h_file(code_param, max_iterations, decoder_type, input_decoder_c, x_hat, detected_data);
          #}
        end

        detected_data = detected_data(1:code_param.data_bits_per_frame);
      end

      if strcmp(code, 'golay')
        for f=1:Ncodewordsperframe
          st = (f-1)*code_param.code_bits_per_frame+1; en=st+code_param.code_bits_per_frame-1;
          arx_codeword = egolaydec(real(r(st:en)) < 0);
          rx_bits = [rx_bits arx_codeword(code_param.data_bits_per_frame+1:code_param.code_bits_per_frame)];
        end
      end

      if strcmp(code, 'diversity')
        for f=1:Ncodewordsperframe
          st = (f-1)*Nc+1;
          r_combined = r(st:st+Nc/2-1) + r(st+Nc/2:st+Nc-1);
          arx_data = real(r_combined) < 0;
          rx_bits = [rx_bits arx_data];
        end

        #{
        % simulate low rate code to mop up errors

        error_positions = xor(rx_bits, tx_bits);
        Nerrs = sum(error_positions);
        if Nerrs < 6
          rx_bits = tx_bits;
        end        
        #}
      end

      rx_bits_log = [rx_bits_log rx_bits];

      error_positions = xor(rx_bits, tx_bits);
      Nerrs = sum(error_positions);
      Nerrs_log = [Nerrs_log Nerrs];
      
      if Nerrs>0,  Ferrs += 1;  end
      Terrs = Terrs + Nerrs;
      Tbits = Tbits + code_param.data_bits_per_frame*Ncodewordsperframe;
    end
      
    % Alternative Codec 2 packet rate measurement indep of framesize

    Nerrs_codec2_log = [];
    Ncodecpacketsize = 28;
    Perrs = 0; Npackets = floor(length(tx_bits_log)/Ncodecpacketsize);
    for p=1:Ncodecpacketsize:Npackets*Ncodecpacketsize
      Nerrs = sum(xor(tx_bits_log(p:p+Ncodecpacketsize-1), rx_bits_log(p:p+Ncodecpacketsize-1)));
      if Nerrs
        Perrs++;
      end
      Nerrs_codec2_log = [Nerrs_codec2_log Nerrs];
    end

    printf("Coded EbNo: %3.1f dB BER: %5.4f PER: %5.4f Nbits: %4d Nerrs: %4d Tpackets: %4d Perr: %4d\n", 
            EbNodB, Terrs/Tbits, Ferrs/ Ntrials, Tbits, Terrs, Ntrials, Ferrs);
    EbNodB_raw = EsNodB - 10*log10(code_param.bits_per_symbol);
    printf("Raw EbNo..: %3.1f dB BER: %5.4f Nbits: %4d Nerrs: %4d\n", EbNodB_raw, 
            Terrs_raw/Tbits_raw, Tbits_raw, Terrs_raw);
    printf("Codec 2 PER: %5.4f Npackets: %d Perrs: %d\n", Perrs/Npackets, Npackets, Perrs);
    
    BERvec(ne) = Terrs/Tbits;
    PERvec(ne) = Perrs/Npackets;

    sim_out.BERvec = BERvec;
    sim_out.PERvec = PERvec;

    error_positions = xor(tx_bits_log, rx_bits_log);
    sim_out.error_positions = error_positions;

    if verbose
      figure(3); clf;
      plot(real(r_log),imag(r_log),'+')
      axis([-2 2 -2 2])
      title('Scatter');

      figure(4); clf;
      subplot(211);
      stem((1:Ntrials)*Tp, Nerrs_raw_log);
      subplot(212);
      stem((1:Npackets)*Tp_codec2, Nerrs_codec2_log);

      if hf_en
        figure(6); clf;

        % limit mesh plot to Np points to plot quickly
      
        Np = 500;
        step = ceil(hf_r/Np);
        mesh(1:Nc, (1:step:hf_r-1)/Rs, abs(hf_model(1:step:hf_r-1,:)))
        title('HF channel amplitude');
        xlabel('Carrier');
        ylabel('Time (s)');
      end

    end
  end
endfunction


function plot_curves(Ntrials=100, hf_en=0)

  if hf_en
    epslabel = 'hf';
  else
    epslabel = 'awgn';
  end

  sim_in.genie_Es    = 1;
  sim_in.code        = 'ldpc';
  sim_in.hf_en       = hf_en;
  sim_in.verbose     = 0;

  if hf_en
    EbNovec = 2:0.5:10; 
  else
    EbNovec = 0:0.5:8; 
  end
  sim_in.EbNovec = EbNovec;

  load HRA_112_112.txt
  load HRA_112_56.txt
  load HRA_56_56.txt
  load HRA_56_28.txt

  printf("HRA_112_112-------------\n");
  sim_out1 = run_sim(sim_in, HRA_112_112, Ntrials);

#{
  printf("HRA_112_56-------------\n");
  sim_out2 = run_sim(sim_in, HRA_112_56 , Ntrials);
#}

  printf("HRA_56_56-------------\n");
  sim_out3 = run_sim(sim_in, HRA_56_56  , Ntrials*2);
#{
  printf("HRA_56_28-------------\n");
  sim_out4 = run_sim(sim_in, HRA_56_28  , Ntrials*2);

  printf("Golay -------------\n");
  sim_in.code = 'golay';
  sim_out5 = run_sim(sim_in, [], Ntrials);
  
#}
  printf("Diversity -------------\n");
  sim_in.code = 'diversity';
  sim_out6 = run_sim(sim_in, [], Ntrials);

  if hf_en
    Ebvec_theory = 2:0.5:12;
    EbNoLin = 10.^(Ebvec_theory/10);
    uncoded_BER_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));
  else
    Ebvec_theory = 0:0.5:8;
    uncoded_BER_theory = 0.5*erfc(sqrt(10.^(Ebvec_theory/10)));
  end

  % need standard packet size to compare
  % packet error if bit 0, or bit 1, or bit 2 .....
  %              or bit 0 and bit 1
  % no packet error if all bits ok (1-p(0))*(1-p(1))
  % P(packet error) = p(0)+p(1)+....

  uncoded_PER_theory = 1 - (1-uncoded_BER_theory).^28;

  figure(1); clf;
  semilogy(Ebvec_theory,  uncoded_BER_theory, 'b+-;BPSK theory;','markersize', 10, 'linewidth', 2)
  hold on;
  semilogy(EbNovec, sim_out1.BERvec, 'g+-;rate 1/2 HRA 112 112;','markersize', 10, 'linewidth', 2)
  %semilogy(EbNovec, sim_out2.BERvec, 'r+-;rate 2/3 HRA 112 56;','markersize', 10, 'linewidth', 2)
  semilogy(EbNovec, sim_out3.BERvec, 'c+-;rate 1/2 HRA 56 56;','markersize', 10, 'linewidth', 2)
  %semilogy(EbNovec, sim_out4.BERvec, 'k+-;rate 2/3 HRA 56 28;','markersize', 10, 'linewidth', 2)
  %semilogy(EbNovec, sim_out5.BERvec, 'm+-;rate 1/2 Golay (24,12);','markersize', 10, 'linewidth', 2)
  semilogy(EbNovec, sim_out6.BERvec, 'go-;rate 1/2 Diversity;','markersize', 10, 'linewidth', 2)
  hold off;
  xlabel('Eb/No')
  ylabel('BER')
  grid
  legend("boxoff");
  axis([min(Ebvec_theory) max(Ebvec_theory) 1E-3 2E-1]);
  epsname = sprintf("ldpc_short_%s_ber.eps", epslabel);
  print('-deps', '-color', epsname)

  figure(2); clf;
  semilogy(Ebvec_theory,  uncoded_PER_theory, 'b+-;BPSK theory;','markersize', 10, 'linewidth', 2)
  hold on;
  semilogy(EbNovec, sim_out1.PERvec, 'g+-;rate 1/2 HRA 112 112;','markersize', 10, 'linewidth', 2)
  %semilogy(EbNovec, sim_out2.PERvec, 'r+-;rate 2/3 HRA 112 56;','markersize', 10, 'linewidth', 2)
  semilogy(EbNovec, sim_out3.PERvec, 'c+-;rate 1/2 HRA 56 56;','markersize', 10, 'linewidth', 2)
  %semilogy(EbNovec, sim_out4.PERvec, 'k+-;rate 2/3 HRA 56 28;','markersize', 10, 'linewidth', 2)
  %semilogy(EbNovec, sim_out5.PERvec, 'm+-;rate 1/2 Golay (24,12);','markersize', 10, 'linewidth', 2)
  semilogy(EbNovec, sim_out6.PERvec, 'go-;rate 1/2 Diversity;','markersize', 10, 'linewidth', 2)
  hold off;
  xlabel('Eb/No')
  ylabel('PER')
  grid
  legend("boxoff");
  if hf_en
    legend("location", "southwest");
  end
  axis([min(Ebvec_theory) max(Ebvec_theory) 1E-2 1]);
  epsname = sprintf("ldpc_short_%s_per.eps", epslabel);
  print('-deps', '-color', epsname)
endfunction


function run_single(bits, code = 'ldpc', channel = 'awgn', EbNodB, interleave=0, error_pattern_filename)

  sim_in.code = code;
  load HRA_112_112.txt
  HRA = HRA_112_112;

  sim_in.genie_Es = 1;

  sim_in.EbNovec = EbNodB;
  if strcmp(channel, 'awgn')
    sim_in.hf_en = 0;
  else
    sim_in.hf_en = 1;
  end

  if interleave
    sim_in.interleave_en = 1;
    sim_in.interleave_frames = interleave;
  else
    sim_in.interleave_frames = 1;
  end
  Ntrials = floor(bits/(112*sim_in.interleave_frames));
  sim_in.verbose = 1;
  sim_out = run_sim(sim_in, HRA, Ntrials);

  if nargin == 6
    fep = fopen(error_pattern_filename, "wb");
    fwrite(fep, sim_out.error_positions, "short");
    fclose(fep);
  end

endfunction


% Used to generate C header file for C port

function run_c_header

  sim_in.code = 'ldpc';
  load HRA_112_112.txt
  data_bits_per_frame = 112;
  rate = 0.5;
  bits = data_bits_per_frame;
  Ntrials = bits/data_bits_per_frame;
  sim_in.genie_Es    = 1;
  sim_in.hf_en = 0;
  sim_in.Esvec = 2;
  sim_in.c_include_file = "../src/HRA_112_112.h";

  sim_out = run_sim(sim_in, HRA_112_112, Ntrials);
endfunction


% Start simulation here ----------------------------------------------

% change this path for your CML installation

init_cml('/home/david/Desktop/cml');

more off;
format;

% simple single point test

run_single(700*60, 'ldpc', 'hf', 6, 32)

% plotting curves (may take a while)

%plot_curves(0);
%plot_curves(500,1);

% generating error files 

%run_single(700*10, 'ldpc', 'awgn', 3, 0, 'awgn_3dB_ldpc.err')
%run_single(700*10, 'diversity', 'awgn', 3, 0, 'awgn_3dB_diversity.err')
%run_single(700*10, 'ldpc', 'hf', 10, 0, 'hf_10dB_ldpc.err')
%run_single(700*10, 'diversity', 'hf', 10, 0, 'hf_10dB_diversity.err')

% generate C header for C port of code

%run_c_header
