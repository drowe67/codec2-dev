% tqpsk.m
%
% David Rowe Sep 2015
%
% QPSK modem Octave test scripts

% TODO:
%   correct freq offset
%   output tx frame to a file of bits for 3rd party modulator
%   improve phase est impl loss using longer window
%   separate into mod and demod functions
%   generate filter coeffs, save to C header file
%   fine timing estimator and meas IL
%   coarse timing (frame sync)

% make sure we get same random numbers every time

more off;
qpsk;
ldpc;


function [aqpsk code_param] = init_modem_ldpc

  % Init modem ------------------------------------------------------------------------------

  aqpsk.Rbdata        = 2400;                                % bit rate (Hz) of payload data
  aqpsk.code_rate     = 0.5;
  aqpsk.framesize     = 1152;                                % number of bits in LDPC encoded frame (data plus parity bits)
  aqpsk.Ndatabits     = aqpsk.framesize*aqpsk.code_rate;
  aqpsk.Rbcoded       = aqpsk.Rbdata/aqpsk.code_rate;        % bit rate (Hz) of LPDC coded data
  aqpsk.Rf            = aqpsk.Rbcoded/aqpsk.framesize;       % frame rate (Hz)
  aqpsk.bits_per_symb = 2;                                   % QPSK has two bits per symbol
  aqpsk.Rscoded       = aqpsk.Rbcoded/aqpsk.bits_per_symb;   % symbol rate of LDPC coded data
  aqpsk.Nsymb         = aqpsk.framesize/aqpsk.bits_per_symb; % number of coded data symbols/frame
  aqpsk.Npilotstep    = 8;                                   % number of data symbols between each pilot symbol

  aqpsk.Npilotsframe  = (aqpsk.framesize/aqpsk.bits_per_symb)/aqpsk.Npilotstep;
  if aqpsk.Npilotsframe != floor(aqpsk.Npilotsframe)
    error("Npilotsframe must be an integer");
  end

  aqpsk.Nsymbpilot    = aqpsk.Nsymb + aqpsk.Npilotsframe;   % total number of symbols per frame including pilots
  aqpsk.Rs            = aqpsk.Rf*aqpsk.Nsymbpilot;          % total symbol rate (Hz)
  aqpsk.pilots        = ones(1, aqpsk.Npilotsframe);
  aqpsk.Np            = 2;                                  % number of pilots used for phase estimation
  aqpsk.rx_pilot_buf  = zeros(1,aqpsk.Np);

  % Init LDPC --------------------------------------------------------------------------

  currentdir = pwd;
  addpath '/home/david/tmp/cml/mat'    % assume the source files stored here
  cd /home/david/tmp/cml
  CmlStartup                           % note that this is not in the cml path!
  cd(currentdir)

  mod_order = 4; 
  modulation = 'QPSK';
  mapping = 'gray';
  demod_type = 0;
  decoder_type = 0;
  max_iterations = 100;

  code_param = ldpc_init(aqpsk.code_rate, aqpsk.framesize, modulation, mod_order, mapping);
  assert(code_param.data_bits_per_frame == aqpsk.Ndatabits);

endfunction


function sim_out = run_simulation(sim_in)
  Ntrials      = sim_in.Ntrials;
  EsNodBvec    = sim_in.EsNodBvec;
  verbose      = sim_in.verbose;
  phase_offset = sim_in.phase_offset;
  phase_est_en = sim_in.phase_est_en;
  verbose      = sim_in.verbose;

  rand('state',1); 
  randn('state',1);

  % Init modem and LPDC ----------------------------------------------------------------

  [aqpsk code_param] = init_modem_ldpc;

  Nsymb = aqpsk.Nsymb;
  Nsymbpilot = aqpsk.Nsymbpilot;

  % Start simulation -------------------------------------------------------------------

  for ne = 1:length(EsNodBvec)
    EsNodB = EsNodBvec(ne);
    EsNo = 10^(EsNodB/10);
    variance = 1/EsNo;

    Tbits = Terrs_uc = Terrs = Ferrs = 0;

    % Encode a bunch of frames

    for nn = 1: Ntrials        
      tx_bits = round(rand(1,  aqpsk.Ndatabits));
      [tx_codeword, tx_coded_symb] = ldpc_enc(tx_bits, code_param);
      symbpilot_tx = insert_pilots(tx_coded_symb, aqpsk.pilots, aqpsk.Npilotstep);

      noise = sqrt(variance*0.5)*(randn(1,aqpsk.Nsymbpilot) + j*randn(1,aqpsk.Nsymbpilot));
      symbpilot_rx = symbpilot_tx * exp(j*phase_offset) + noise;

      if phase_est_en
        symbpilot_rx = correct_phase_offset(aqpsk, symbpilot_rx);
      end
      rx_coded_symb = remove_pilots(symbpilot_rx, aqpsk.pilots, aqpsk.Npilotstep);
      rx_codeword = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, rx_coded_symb, EsNo);

      % measure coded and uncoded errors

      errors_positions = xor(tx_bits, rx_codeword(1:aqpsk.Ndatabits));
      Nerrs = sum(errors_positions);

      rx_bits_uncoded = zeros(1, aqpsk.Ndatabits);
      for s=1:Nsymb*aqpsk.code_rate;
         rx_bits_uncoded(2*(s-1)+1:2*s) = qpsk_demod(rx_coded_symb(s));
      end
      errors_positions = xor(tx_bits, rx_bits_uncoded);
      Nerrs_uc = sum(errors_positions);

      if verbose == 2
        % print "." if frame decoded without errors, 'x' if we can't decode

        if Nerrs > 0, printf('x'),  else printf('.'),  end
        if (rem(nn, 50)==0),  printf('\n'),  end    
      end

      if Nerrs > 0,  Ferrs = Ferrs + 1;  end
      Terrs = Terrs + Nerrs;
      Terrs_uc = Terrs_uc + Nerrs_uc;
      Tbits = Tbits + aqpsk.Ndatabits;        
    end

    if verbose
      printf("EsNodB: %3.2f BER: %f Tbits: %d Terrs: %d Ferrs: %d Terrs_uc: %d\n", EsNodB, Terrs/Tbits, Tbits, Terrs,  Ferrs, Terrs_uc);
    end

    sim_out.Tbits(ne)     = Tbits;
    sim_out.Terrs(ne)     = Terrs;
    sim_out.Terrs_uc(ne)  = Terrs_uc;
    sim_out.Ferrs(ne)     = Ferrs;
    sim_out.BER(ne)       = Terrs/Tbits;
    sim_out.FER(ne)       = Ferrs/Ntrials;
    sim_out.BER_uc(ne)    = Terrs_uc/Tbits;
  end

  sim_out.code_rate = aqpsk.code_rate;
end


% ----------------------------------------------------------------------


function symbol_rate_curves
  sim_in.phase_est_en = 1;
  sim_in.phase_offset = 0;
  sim_in.verbose      = 1;
  sim_in.EsNodBvec    = -2:10;
  sim_in.Ntrials      = 10;

  sim_out = run_simulation(sim_in);

  EbNodB = sim_in.EsNodBvec - 10*log10(2);  % QPSK has two bits/symbol
  uncoded_BER_theory = 0.5*erfc(sqrt(10.^(EbNodB/10)));

  figure(1)
  clf
  semilogy(EbNodB, uncoded_BER_theory,'r;uncoded QPSK theory;')
  hold on;
  semilogy(EbNodB, sim_out.BER_uc,'r+;uncoded QPSK simulated;')
  semilogy(EbNodB-10*log10(sim_out.code_rate), sim_out.BER+1E-10,'g;LDPC coded QPSK simulated;');
  hold off;
  grid('minor')
  xlabel('Eb/No (dB)')
  ylabel('BER')
  axis([min(EbNodB) max(EbNodB) min(uncoded_BER_theory) 1])

endfunction


% We are using a 3rd party modulator, this function generates bits
% that can be passed to it.  It's organised as a file of 8 bit bytes
% MSB is sent first.  Each pair of bits is mapped to phases as:
%
% MSB LSB   Phase
%   0   0       0
%   0   1      90
%   1   1     180
%   1   0     270 

function generate_modulator_input_file(filename, seconds)
         
  [aqpsk code_param] = init_modem_ldpc;

  Nframes = floor(seconds/aqpsk.Rf);
  Nsecs = Nframes*aqpsk.Rf;
  Nbitsframe = aqpsk.Nsymbpilot*aqpsk.bits_per_symb;
  Nbytesframe = Nbitsframe/8;

  printf("Frames: %d Seconds: %3.1f QPSK Symb Rate: %5.2f\n", Nframes, Nsecs, aqpsk.Rs);
  printf("Bit/frame: %d  Packed Bytes/frame: %4.1f\n", Nbitsframe, Nbytesframe);

  fbit = fopen(filename,"wb"); 
 
  for nn = 1:Nframes
    tx_bits = round(rand(1,  aqpsk.Ndatabits));
    [tx_codeword, tx_coded_symb] = ldpc_enc(tx_bits, code_param);
    symbpilot_tx = insert_pilots(tx_coded_symb, aqpsk.pilots, aqpsk.Npilotstep);

    % demodulate QPSK symbols back to bits
   
    tx_frame = zeros(1, aqpsk.Nsymbpilot*aqpsk.bits_per_symb);
    for s=1:aqpsk.Nsymbpilot
      tx_frame(2*(s-1)+1:2*s) = qpsk_demod(symbpilot_tx(s));
    end
    packed_tx_frame = packmsb(tx_frame);
    fwrite(fbit, packed_tx_frame, "uchar");
  end
  
  fclose(fbit);

endfunction

% ---------------------------------------------------------------------------------

generate_modulator_input_file("modinputbits_v0.1.bin", 300)
