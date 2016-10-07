% test_ldpc_fsk_lib
% David Rowe 16 April 2016
%
% Some tests for ldpc_fsk_lib, and C versions ldpc_enc and ldpc_dec

1;

% encodes and decodes one frame, also writes codeword.bin for testing
% decode_from_file() below, and can optionally generate include file for
% C version of decoder.

function [data code_param] = simple_ut(c_include_file)
  load('H2064_516_sparse.mat');
  HRA = full(HRA);  
  max_iterations = 100;
  decoder_type = 0;
  EsNodB = 3;
  mod_order = 2;

  code_param = ldpc_init(HRA, mod_order);
  data = round( rand( 1, code_param.data_bits_per_frame ) );
  codeword = ldpc_encode(code_param, data);
  f = fopen("codeword.bin","wt"); fwrite(f, codeword, "uint8"); fclose(f);
  s = 1 - 2 * codeword;   
  code_param.symbols_per_frame = length( s );

  EsNo = 10^(EsNodB/10);
  variance = 1/(2*EsNo);
  noise = sqrt(variance)* randn(1,code_param.symbols_per_frame); 
  rx = s + noise;
  
  if nargin == 1
    code_param.c_include_file = c_include_file;
  end
  [detected_data Niters] = ldpc_decode(rx, code_param, max_iterations, decoder_type);
  
  error_positions = xor(detected_data(1:code_param.data_bits_per_frame), data);
  Nerrs = sum(error_positions);

  printf("Nerrs = %d\n", Nerrs);
end


% This version decodes from a file of bits

function detected_data = decode_from_file(filename)
  max_iterations = 100;
  decoder_type = 0;
  load('H2064_516_sparse.mat');
  HRA = full(HRA);  
  mod_order = 2;

  f = fopen(filename,"rb"); codeword = fread(f, "uint8")'; fclose(f);
  r = 1 - 2 * codeword;   
  code_param = ldpc_init(HRA, mod_order);
  [detected_data Niters] = ldpc_decode(r, code_param, max_iterations, decoder_type);
end


% plots a BER curve for the decoder.  Takes a while to run, uses parallel cores

function plot_curve
  num_cores = 4;              % set this to the number of cores you have

  load('H2064_516_sparse.mat');
  HRA = full(HRA);  
  [Nr Nc] = size(HRA); 
  sim_in.rate = (Nc-Nr)/Nc;

  sim_in.HRA            = HRA;
  sim_in.mod_order      = 2;
  sim_in.framesize      = Nc;
  sim_in.mod_order      = 2; 
  sim_in.Lim_Ferrs      = 100;

  % note we increase number of trials as BER goes down

  Esvec   = [   0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 ]; 
  Ntrials = [ 1E4 1E4 1E4 1E4 1E5 1E5 1E5 1E5 1E5 ];
  num_runs = length(Esvec)

  sim_in_vec(1:num_runs) = sim_in;
  for i = 1:num_runs
    sim_in_vec(i).Esvec   = Esvec(i);
    sim_in_vec(i).Ntrials = Ntrials(i);
  end

  %sim_out = ldpc5(sim_in_vec(1));
  tstart = time();
  sim_out = pararrayfun(num_cores, @ldpc5, sim_in_vec);
  tend = time();

  total_bits = sum(Ntrials)*sim_in.framesize;
  total_secs = tend - tstart;
  printf("%d bits in %4.1f secs, or %5f bits/s\n", total_bits, total_secs, total_bits/total_secs);

  for i=1:num_runs
    Ebvec(i)  = sim_out(i).Ebvec;
    BERvec(i) = sim_out(i).BERvec;
  end
  semilogy(Ebvec,  BERvec, '+-')
  xlabel('Eb/N0')
  ylabel('BER')
  title(['H2064 516 sparse.mat' ' ' date])

end


% Test C encoder

function test_c_encoder
  load('H2064_516_sparse.mat');
  HRA = full(HRA);  
  max_iterations = 100;
  decoder_type = 0;
  EsNodB = 3;
  mod_order = 2;
  frames = 100;

  EsNo = 10^(EsNodB/10);
  variance = 1/(2*EsNo);

  code_param = ldpc_init(HRA, mod_order);

  data = round(rand(1,frames*code_param.data_bits_per_frame));
  f = fopen("data.bin","wt"); fwrite(f, data, "uint8"); fclose(f);

  % Outboard C encoder

  system("../src/ldpc_enc data.bin codewords.bin"); 

  % Test with Octave decoder

  f = fopen("codewords.bin","rb"); codewords = fread(f, "uint8")'; fclose(f);
  
  Nerrs = 0;
  for i=1:frames
    st = (i-1)*code_param.symbols_per_frame+1; en = st+code_param.symbols_per_frame-1;
    tx = 1 - 2 * codewords(st:en);   

    noise = sqrt(variance)*randn(1,code_param.symbols_per_frame); 
    rx = tx + noise;

    [detected_data Niters] = ldpc_decode(rx, code_param, max_iterations, decoder_type);

    st = (i-1)*code_param.data_bits_per_frame+1; en = st+code_param.data_bits_per_frame-1;
    error_positions = xor(detected_data(1:code_param.data_bits_per_frame), data(st:en));
    Nerrs += sum(error_positions);
  end

  printf("Nerrs = %d\n", Nerrs);
end


function test_c_decoder
  load('H2064_516_sparse.mat');
  HRA = full(HRA);  
  max_iterations = 100;
  decoder_type = 0;
  mod_order = 2;
  frames = 10;
  EsNodB = 2;
  sdinput = 1;

  EsNo = 10^(EsNodB/10);
  variance = 1/(2*EsNo);

  code_param = ldpc_init(HRA, mod_order);
  data = round(rand(1,code_param.data_bits_per_frame*frames));
  
  f = fopen("data.bin","wt"); fwrite(f, data, "uint8"); fclose(f);
  system("../src/ldpc_enc data.bin codewords.bin"); 
  f = fopen("codewords.bin","rb"); codewords = fread(f, "uint8")'; fclose(f);

  s = 1 - 2 * codewords;   
  noise = sqrt(variance)*randn(1,code_param.symbols_per_frame*frames); 
  r = s + noise;

  % calc LLRs frame by frame

  for i=1:frames
    st = (i-1)*code_param.symbols_per_frame+1;
    en = st + code_param.symbols_per_frame-1;
    llr(st:en) = sd_to_llr(r(st:en));    
  end

  % Outboard C decoder

  if sdinput
    f = fopen("sd.bin","wb"); fwrite(f, r, "double"); fclose(f);
    system("../src/ldpc_dec sd.bin data_out.bin --sdinput"); 
  else
    f = fopen("llr.bin","wb"); fwrite(f, llr, "double"); fclose(f);
    system("../src/ldpc_dec llr.bin data_out.bin"); 
  end

  f = fopen("data_out.bin","rb"); data_out = fread(f, "uint8")'; fclose(f);
  
  Nerrs = Nerrs2 = zeros(1,frames);
  for i=1:frames

    % Check C decoder
    
    data_st = (i-1)*code_param.data_bits_per_frame+1;
    data_en = data_st+code_param.data_bits_per_frame-1;
    st = (i-1)*code_param.symbols_per_frame+1;
    en = st+code_param.data_bits_per_frame-1;
    data_out_c = data_out(st:en);
    error_positions = xor(data_out_c, data(data_st:data_en));
    Nerrs(i) = sum(error_positions);

    % Octave decoder 

    st = (i-1)*code_param.symbols_per_frame+1; en = st+code_param.symbols_per_frame-1;
    [detected_data Niters] = ldpc_decode(r(st:en), code_param, max_iterations, decoder_type);
    st = (i-1)*code_param.data_bits_per_frame+1; en = st+code_param.data_bits_per_frame-1;
    data_out_octave = detected_data(1:code_param.data_bits_per_frame);
    error_positions = xor(data_out_octave, data(st:en));
    Nerrs2(i) = sum(error_positions);
    %printf("%4d ", Niters);
  end
  printf("Errors per frame:\nC.....:");
  for i=1:frames
    printf("%4d ", Nerrs(i));
  end
  printf("\nOctave:");
  for i=1:frames
    printf("%4d ", Nerrs2(i));
  end
  printf("\n");

end

% Saves a complex vector s to a file "filename" of IQ unsigned 8 bit
% chars, same as RTLSDR format.

function save_rtlsdr(filename, s)
  mx = max(abs(s));
  re = real(s); im = imag(s);
  l = length(s);
  iq = zeros(1,2*l);
  iq(1:2:2*l) = 127 + re*(127/mx); 
  iq(2:2:2*l) = 127 + im*(127/mx); 
  fs = fopen(filename,"wb");
  fwrite(fs,iq,"uchar");
  fclose(fs);
endfunction


% Using simulated SSTV packet, generate complex fsk mod signals, 8-bit
% unsigned IQ for feeding into demod chain
#{
todo: [X] uncoded BER
          [X] octave fsk demod
          [X] use C demod
      [ ] compared unsigned 8 bit IQ to regular 16-bit
          + modulate complex signal
          + check BER
          + feed into csdr stack
      [ ] test with resampler
      [ ] measure effect on PER with coding
#}

function [n_uncoded_errs n_uncoded_bits] = run_sstv_sim(sim_in, EbNodB)

  frames = sim_in.frames;
  c_demod = sim_in.c_demod;

  % init LDPC code

  load('H2064_516_sparse.mat');
  HRA = full(HRA);  
  max_iterations = 100;
  decoder_type = 0;
  mod_order = 2;

  code_param = ldpc_init(HRA, mod_order);

  % note fixed frame of bits used for BER testing

  tx_codeword = gen_sstv_frame;

  % init FSK modem

  fsk_horus_as_a_lib = 1;
  fsk_horus;
  states         = fsk_horus_init_hbr(9600, 8, 1200, 2, length(tx_codeword));
  states.ftx     = [1200 2400];
  states.df(1:states.M) = 0;
  states.dA(1:states.M) = 1;
  states.tx_real = 0;

  % set up AWGN channel 

  EbNo = 10^(EbNodB/10);
  variance = states.Fs/(states.Rs*EbNo*states.bitspersymbol);

  % start simulation ----------------------------------------

  tx_bit_stream = [];
  for i=1:frames+1
    tx_bit_stream = [tx_bit_stream tx_codeword];
  end

  % modulate and channel model

  tx = fsk_horus_mod(states, tx_bit_stream);
  noise_real = sqrt(variance)*randn(length(tx),1);
  noise_complex = sqrt(variance/2)*(randn(length(tx),1) + j*randn(length(tx),1));

  % demodulate

  if c_demod
    if states.tx_real
      rx = tx + noise_real;
    else
      rx = 2*real(tx) + noise_real;
    end
    SNRdB = 10*log10(var(tx)/var(noise_real));
    rx_scaled = 1000*real(rx);
    f = fopen("fsk_demod.raw","wb"); fwrite(f, rx_scaled, "short"); fclose(f);
    system("../build_linux/src/fsk_demod 2X 8 9600 1200 fsk_demod.raw fsk_demod.bin");
    f = fopen("fsk_demod.bin","rb"); rx_bit_stream = fread(f, "uint8")'; fclose(f);
  else
    if states.tx_real
      rx = tx + noise_real;
    else
      rx = tx + noise_complex;
    end
    SNRdB = 10*log10(var(tx)/var(noise_complex));

    % demodulate frame by frame using Octave demod

    st = 1;
    run_frames = floor(length(rx)/states.N);
    rx_bit_stream = [];
    rx_sd_stream = [];
    for f=1:run_frames

      % extract nin samples from rx sample stream

      nin = states.nin;
      en = st + states.nin - 1;

      if en <= length(rx) % due to nin variations its possible to overrun buffer
        sf = rx(st:en);
        st += nin;

        % demodulate to stream of bits

        states.f = [1200 2400];
        [rx_bits states] = fsk_horus_demod(states, sf);
        rx_bit_stream = [rx_bit_stream rx_bits];
        rx_sd_stream = [rx_sd_stream states.rx_bits_sd];
      end
    end
  end

  % state machine. Look for SSTV UW.  When found count bit errors over one frame of bits

  state = "wait for uw";
  start_uw_ind = 16*10+1; end_uw_ind = start_uw_ind + 2*10 - 1;
  uw_rs232 = tx_codeword(start_uw_ind:end_uw_ind); luw = length(uw_rs232);
  start_frame_ind =  end_uw_ind + 1;
  nbits = length(rx_bit_stream);
  uw_thresh = 0;
  n_uncoded_errs = 0;
  n_uncoded_bits = 0;

  % might as well include RS232 framing bits in uncoded error count

  nbits_frame = code_param.data_bits_per_frame*10/8;  

  for i=luw:nbits
    next_state = state;
    if strcmp(state, 'wait for uw')
      uw_errs = xor(rx_bit_stream(i-luw+1:i), uw_rs232);
      if uw_errs <= uw_thresh
        next_state = 'count errors';
        tx_frame_ind = start_frame_ind;
        rx_frame_ind = i + 1;
        n_uncoded_errs_this_frame = 0;
        %printf("%d %s %s\n", i, state, next_state);
      end
    end
    if strcmp(state, 'count errors')
      n_uncoded_errs_this_frame += xor(rx_bit_stream(i), tx_codeword(tx_frame_ind));
      n_uncoded_bits++;
      tx_frame_ind++;
      if tx_frame_ind == (start_frame_ind+nbits_frame)
        n_uncoded_errs += n_uncoded_errs_this_frame;
        printf("n_uncoded_errs_this_frame: %d\n", n_uncoded_errs_this_frame);
        frame_rx232_rx = rx_bit_stream(rx_frame_ind:rx_frame_ind+nbits_frame-1);
        %tx_codeword(start_frame_ind+1:start_frame_ind+10)
        %frame_rx232_rx(1:10)
        sstv_checksum(frame_rx232_rx);
        next_state = 'wait for uw';
      end
    end
    state = next_state;
  end

  uncoded_ber = n_uncoded_errs/n_uncoded_bits;
  printf("EbNodB: %4.1f SNRdB: %4.1f n_uncoded_bits: %d n_uncoded_errs: %d BER: %4.3f\n", 
          EbNodB, SNRdB, n_uncoded_bits, n_uncoded_errs, uncoded_ber);  
endfunction

% Start simulation --------------------------------------------------------

more off;
currentdir = pwd;
thiscomp = computer;

if strcmpi(thiscomp, 'MACI64')==1
   if exist('CMLSimulate')==0
        cd '/Users/bill/Current/Projects/DLR_FSO/Visit2013_FSO_GEO/cml'
        addpath '../'    % assume the source files stored here
        CmlStartup       % note that this is not in the cml path!
        disp('added MACI64 path and run CmlStartup')
    end
end

if strfind(thiscomp, 'pc-linux-gnu')==8 
   if exist('LdpcEncode')==0, 
        cd '~/tmp/cml'
        CmlStartup 
        disp('CmlStartup has been run')
	% rmpath '/home/bill/cml/mexhelp'  % why is this needed? 
	% maybe different path order in octave cf matlab ? 
    end
end

cd(currentdir)

ldpc_fsk_lib;
randn('state',1);
rand('state',1);

% ------------------ select which demo/test to run here ---------------

demo = 10;

if demo == 1
  printf("simple_ut....\n");
  data = simple_ut;
end

if demo == 2
  printf("generate C header file....\n");
  data = simple_ut("../src/H2064_516_sparse.h");
end

if demo == 3
  printf("decode_from_file ......\n");
  data = simple_ut;
  detected_data = decode_from_file("codeword.bin");
  error_positions = xor( detected_data(1:length(data)), data );
  Nerrs = sum(error_positions);
  printf("  Nerrs = %d\n", Nerrs);
end

if demo == 4
  printf("plot a curve....\n");
  plot_curve;
end

if demo == 5

  % generate test data and save to disk

  [data code_param] = simple_ut;
  f = fopen("dat_in2064.bin","wb"); fwrite(f, data, "uint8"); fclose(f);

  % Outboard C encoder

  system("../src/ldpc_enc dat_in2064.bin dat_op2064.bin"); 

  % Test with Octave decoder

  detected_data = decode_from_file("dat_op2064.bin");
  error_positions = xor(detected_data(1:length(data)), data);
  Nerrs = sum(error_positions);
  printf("Nerrs = %d\n", Nerrs);
end

if demo == 6
  test_c_encoder;
end

if demo == 7
  test_c_decoder;
end

if demo == 8
  frames = 100;
  EsNodB = 3;
  EsNo = 10^(EsNodB/10);
  variance = 1/(2*EsNo);

  frame_rs232 = [];
  for i=1:frames
    frame_rs232 = [frame_rs232 gen_sstv_frame];
  end

  % write hard decn version to disk file, useful for fsk_mod input

  f = fopen("sstv.bin","wb"); fwrite(f, frame_rs232, "char"); fclose(f);

  % soft decision version (with noise)

  s = 1 - 2*frame_rs232;
  noise = sqrt(variance)*randn(1,length(frame_rs232)); 
  r = s + noise;
  f = fopen("sstv_sd.bin","wb"); fwrite(f, r, "float32"); fclose(f);
end


if demo == 9
  frames = 100;
  EbNodB = 11;

  frame_rs232 = [];
  for i=1:frames
    frame_rs232 = [frame_rs232 gen_sstv_frame];
  end

  % Use C FSK modulator to generate modulated signal

  f = fopen("sstv.bin","wb"); fwrite(f, frame_rs232, "char"); fclose(f);
  system("../build_linux/src/fsk_mod 2 9600 1200 1200 2400 sstv.bin fsk_mod.raw");

  % Add some channel noise here in Octave

  f = fopen("fsk_mod.raw","rb"); tx = fread(f, "short")'; fclose(f); tx_pwr = var(tx);
  Fs = 9600; Rs=1200; EbNolin = 10 ^ (EbNodB/10);
  variance = (tx_pwr/2)*states.Fs/(states.Rs*EbNolin*states.bitspersymbol);
  noise = sqrt(variance)*randn(1,length(tx)); 
  SNRdB = 10*log10(var(tx)/var(noise));
  rx = tx + noise;
  f = fopen("fsk_demod.raw","wb"); tx = fwrite(f, rx, "short"); fclose(f);
 
  % Demodulate using C modem and C de-framer/LDPC decoder

  system("../build_linux/src/fsk_demod 2XS 8 9600 1200 fsk_demod.raw - | ../src/drs232_ldpc - dummy_out.bin");
end


if demo == 10
  sim_in.frames = 3;
  EbNodBvec = 9;
  
  sim_in.c_demod = 0;
  ber_octave = [];
  for i = 1:length(EbNodBvec)
    [n_uncoded_errs n_uncoded_bits] = run_sstv_sim(sim_in, EbNodBvec(i));
    ber_octave(i) = n_uncoded_errs/n_uncoded_bits;
  end

  sim_in.c_demod = 1;
  ber_c = [];
  for i = 1:length(EbNodBvec)
    [n_uncoded_errs n_uncoded_bits] = run_sstv_sim(sim_in, EbNodBvec(i));
    ber_c(i) = n_uncoded_errs/n_uncoded_bits;
  end

  figure(1);
  clf;
  semilogy(EbNodBvec,  ber_octave, '+-;Octave demod;')
  grid;
  xlabel('Eb/No (dB)')
  ylabel('BER')

  hold on;
  semilogy(EbNodBvec,  ber_c, 'g+-;C demod;')
  legend("boxoff");
  hold off;
end

