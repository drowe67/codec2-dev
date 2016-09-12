% test_ldpc_fsk_lib
% David Rowe 16 April 2016
%
% Some tests for ldpc_fsk_lib

1;

% encodes and decodes one frame, also writes codeword.bin for testing
% decode_from_file() below, and can optionally generate include file for
% C version of decoder.

function [data code_param] = simple_ut(c_include_file)
  load('H2064_516_sparse.mat');
  HRA = full(HRA);  
  max_iterations = 100;
  decoder_type = 0;
  EsNodB = 4;
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


% plots a BER cuvre for the decoder.  Takes a while to run, uses parallel cores

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

% binary flags to run various demos, e.g. "15" to run 1 .. 8

demos = 32;

if bitand(demos,1)
  printf("simple_ut....\n");
  data = simple_ut;
end

if bitand(demos,2)
  printf("generate C header file....\n");
  data = simple_ut("../src/ldpc_code.h");
end

if bitand(demos,4)
  printf("decode_from_file ......\n");
  data = simple_ut;
  detected_data = decode_from_file("codeword.bin");
  error_positions = xor( detected_data(1:length(data)), data );
  Nerrs = sum(error_positions);
  printf("  Nerrs = %d\n", Nerrs);
end

if bitand(demos,8)
  printf("plot a curve....\n");
  plot_curve;
end

if bitand(demos,16)

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


if bitand(demos,32)
  test_c_encoder;
end
