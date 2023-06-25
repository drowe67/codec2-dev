% ofdm_helper.m
%
% Misc functions that are used to support OFDM modem development, that
% aren't required for modem operation

1;

%------------------------------------------------------------------------------
% print_config - utility function to use ascii-art to describe the modem frame
%------------------------------------------------------------------------------

function print_config(states)
  ofdm_load_const;

  % ASCII-art packet visualisation
  s=1; u=1; Nuwsyms=length(uw_ind_sym);
  cr = 1:Nc+2;
  for f=1:Np
    for r=1:Ns
      for c=cr
        if r == 1
          if (c==1) && states.edge_pilots
            sym="P";
          elseif (c==Nc+1) && states.edge_pilots
            sym="P";
          elseif c>1 && c <=(Nc+1)
            sym="P";
          else
            sym=" ";
          end
        elseif c>1 && c <=(Nc+1)
          sym=".";
          if (u <= Nuwsyms) && (s == uw_ind_sym(u)) sym="U"; u++; end
          s++;
        else
          sym=" ";
        end
        printf("%s",sym);
      end
      printf("\n");
    end
  end

  printf("Nc=%d Ts=%4.3f Tcp=%4.3f Ns: %d Np: %d\n", Nc, 1/Rs, Tcp, Ns, Np);
  printf("Nsymperframe: %d Nbitsperpacket: %d Nsamperframe: %d Ntxtbits: %d Nuwbits: %d Nuwframes: %d\n",
          Ns*Nc, Nbitsperpacket, Nsamperframe, Ntxtbits, Nuwbits, Nuwframes);
  printf("uncoded bits/s: %4.1f\n",  Nbitsperpacket*Fs/(Np*Nsamperframe));
end

%-----------------------------------------------------------------------
% create_ldpc_test_frame - generate a test frame of bits
%-----------------------------------------------------------------------

function [tx_bits payload_data_bits codeword] = create_ldpc_test_frame(states, coded_frame=1)
  ofdm_load_const;
  ldpc;
  gp_interleaver;

  if coded_frame
    % Set up LDPC code

    mod_order = 4; bps = 2; modulation = 'QPSK'; mapping = 'gray';

    init_cml(); % TODO: make this path sensible and portable
    load HRA_112_112.txt
    [code_param framesize rate] = ldpc_init_user(HRA_112_112, modulation, mod_order, mapping);
    assert(Nbitsperframe == (code_param.coded_bits_per_frame + Nuwbits + Ntxtbits));

    payload_data_bits = round(ofdm_rand(code_param.data_bits_per_frame)/32767);
    codeword = LdpcEncode(payload_data_bits, code_param.H_rows, code_param.P_matrix);
    Nsymbolsperframe = length(codeword)/bps;

    % need all these steps to get actual raw codeword bits at demod ..

    tx_symbols = [];
    for s=1:Nsymbolsperframe
      tx_symbols = [tx_symbols qpsk_mod( codeword(2*(s-1)+1:2*s) )];
    end

    tx_symbols = gp_interleave(tx_symbols);

    codeword_raw = [];
    for s=1:Nsymbolsperframe
      codeword_raw = [codeword_raw qpsk_demod(tx_symbols(s))];
    end
  else
    codeword_raw = round(ofdm_rand(Nbitsperpacket-(Nuwbits+Ntxtbits))/32767);
  end

  % insert UW and txt bits

  tx_bits = assemble_modem_packet(states, codeword_raw, zeros(1,Ntxtbits));
  assert(Nbitsperpacket == length(tx_bits));

endfunction

% automated test

function test_assemble_disassemble(states)
  ofdm_load_const;

  Nsymsperpacket = Nbitsperpacket/bps;
  Ndatabitsperpacket = Nbitsperpacket-(Nuwbits+Ntxtbits);
  Ndatasymsperpacket = Ndatabitsperpacket/bps;
  codeword_bits = round(ofdm_rand(Ndatabitsperpacket)/32767);
  tx_bits = assemble_modem_packet(states, codeword_bits, zeros(1,Ntxtbits));

  tx_syms = zeros(1,Nsymsperpacket);
  for s=1:Nsymsperpacket
    if bps == 2
      tx_syms(s) = qpsk_mod(tx_bits(bps*(s-1)+1:bps*s));
    elseif bps == 4
      tx_syms(s) = qam16_mod(states.qam16,tx_bits(bps*(s-1)+1:bps*s));
    end
  end
  codeword_syms = zeros(1,Ndatasymsperpacket);
  for s=1:Ndatasymsperpacket
    if bps == 2
      codeword_syms(s) = qpsk_mod(codeword_bits(bps*(s-1)+1:bps*s));
    elseif bps == 4
      codeword_syms(s) = qam16_mod(states.qam16,codeword_bits(bps*(s-1)+1:bps*s));
    end
  end

  [rx_uw rx_codeword_syms payload_amps txt_bits] = disassemble_modem_packet(states, tx_syms, ones(1,Nsymsperpacket));
  assert(rx_uw == states.tx_uw);
  Ndatasymsperframe = (Nbitsperpacket-(Nuwbits+Ntxtbits))/bps;
  assert(codeword_syms == rx_codeword_syms);
endfunction

% test function, kind of like a CRC for QPSK symbols, to compare two vectors

function acc = test_acc(v)
  sre = 0; sim = 0;
  for i=1:length(v)
    x = v(i);
    re = round(real(x)); im = round(imag(x));
    sre += re; sim += im;
    %printf("%d %10f %10f %10f %10f\n", i, re, im, sre, sim);
  end
  acc = sre + j*sim;
end


% Save test bits frame to a text file in the form of a C array
%
% usage:
%   ofdm_lib; test_bits_ofdm_file
%

function test_bits_ofdm_file
  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 17; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  [test_bits_ofdm payload_data_bits codeword] = create_ldpc_test_frame(states);
  printf("%d test bits\n", length(test_bits_ofdm));

  f=fopen("../src/test_bits_ofdm.h","wt");
  fprintf(f,"/* Generated by test_bits_ofdm_file() Octave function */\n\n");
  fprintf(f,"const int test_bits_ofdm[]={\n");
  for m=1:length(test_bits_ofdm)-1
    fprintf(f,"  %d,\n",test_bits_ofdm(m));
  endfor
  fprintf(f,"  %d\n};\n",test_bits_ofdm(end));

  fprintf(f,"\nconst int payload_data_bits[]={\n");
  for m=1:length(payload_data_bits)-1
    fprintf(f,"  %d,\n",payload_data_bits(m));
  endfor
  fprintf(f,"  %d\n};\n",payload_data_bits(end));

  fprintf(f,"\nconst int test_codeword[]={\n");
  for m=1:length(codeword)-1
    fprintf(f,"  %d,\n",codeword(m));
  endfor
  fprintf(f,"  %d\n};\n",codeword(end));

  fclose(f);

endfunction


% Get rid of nasty unfiltered stuff either side of OFDM signal
% This may need to be tweaked, or better yet made a function of Nc, if Nc changes
%
% usage:
%  ofdm_lib; make_ofdm_bpf(1);

function bpf_coeff = make_ofdm_bpf(write_c_header_file)
  filt_n = 100;
  Fs = 8000;

  bpf_coeff  = fir2(filt_n,[0 900 1000 2000 2100 4000]/(Fs/2),[0.001 0.001 1 1 0.001 0.001]);

  if write_c_header_file
    figure(1)
    clf;
    h = freqz(bpf_coeff,1,Fs/2);
    plot(20*log10(abs(h)))
    grid minor

    % save coeffs to a C header file

    f=fopen("../src/ofdm_bpf_coeff.h","wt");
    fprintf(f,"/* 1000 - 2000 Hz FIR filter coeffs */\n");
    fprintf(f,"/* Generated by make_ofdm_bpf() in ofdm_lib.m */\n");

    fprintf(f,"\n#define OFDM_BPF_N %d\n\n", filt_n);

    fprintf(f,"float ofdm_bpf_coeff[]={\n");
    for r=1:filt_n
      if r < filt_n
        fprintf(f, "  %f,\n",  bpf_coeff(r));
      else
        fprintf(f, "  %f\n};", bpf_coeff(r));
      end
    end
    fclose(f);
  end

endfunction

% Helper function to help design UW error thresholds, in particular for raw
% data modes.  See also https://www.rowetel.com/wordpress/?p=7467
function ofdm_determine_bad_uw_errors(Nuw)
   figure(1); clf;
   
   % Ideally the 10% and 50% BER curves are a long way apart
   
   plot(0:Nuw, binocdf(0:Nuw,Nuw,0.1),';BER=0.1;'); hold on; 
   plot(binocdf(0:Nuw,Nuw,0.5),';BER=0.5;'); 
   
   % Suggested threshold for raw data modes is the 5% probability
   % level for the 50% BER curve.  The pre/post-amble has a low chance
   % of failure.  If it does make an error, then we will have random
   % bits presented as the UW (50% BER in UW). This threshold means
   % there is only a 5% case of random bits being accepted as a valid UW
 
   bad_uw_errors = max(find(binocdf(0:Nuw,Nuw,0.5) <= 0.05))+1; 
   plot([bad_uw_errors bad_uw_errors],[0 1],';bad uw errors;'); hold off; grid
   
   xlabel('bits');
   printf("for Nuw = %d, suggest bad_uw_errors = %d\n", Nuw, bad_uw_errors);
end

% Returns level threshold such that threshold_cdf of the tx magnitudes are 
% beneath that level.  Helper function that can be used to design 
% the clipper level.  See also https://www.rowetel.com/?p=7596
function threshold_level = ofdm_determine_clip_threshold(tx, threshold_cdf)
  Nsteps = 25;
  mx = max(abs(tx));
  cdf = empirical_cdf(mx*(1:Nsteps)/Nsteps,abs(tx));
  threshold_level = find(cdf >= threshold_cdf)(1)*mx/25;
  printf("threshold_cdf: %f threshold_level: %f\n", threshold_cdf, threshold_level);
  figure(1); clf; [hh nn] = hist(abs(tx),Nsteps,1);
  plotyy(nn,hh,mx*(1:Nsteps)/Nsteps,cdf); title('PDF and CDF Estimates'); grid;
end


%  helper function that adds channel simulation and ensures we don't saturate int16 output samples  
function [rx_real rx] = ofdm_channel(states, tx, SNR3kdB, channel, freq_offset_Hz)
  [rx_real rx sigma] = channel_simulate(states.Fs, SNR3kdB, freq_offset_Hz, channel, tx, states.verbose);
    
  % multipath models can lead to clipping of int16 samples
  num_clipped = length(find(abs(rx_real>32767)));
  while num_clipped/length(rx_real) > 0.001
    rx_real /= 2;
    num_clipped = length(find(abs(rx_real>32767)));
    printf("WARNING: output samples clipped, reducing level\n")
  end
endfunction


