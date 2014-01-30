% ldpcenc.m
% David Rowe 20 Dec 2013
% 
% LDPC encoder function. Takes a random data pattern, LDPC Encodes and
% inserts Unique Word (UW) sync bits and ouputs this as a packed
% binary file suitable for the Nc=18 carrier FDMDV modulator,
% fdmdv_mod.  Also produces a "modulated" output file of QPSK
% symbols, suitable for feeding into ldpcdec for testing.

function ldpcenc(filename)

  % Start CML library

  currentdir = pwd;
  addpath '/home/david/tmp/cml/mat'    % assume the source files stored here
  cd /home/david/tmp/cml
  CmlStartup                           % note that this is not in the cml path!
  cd(currentdir)
  
  % Our LDPC library

  ldpc;

  % Start simulation

  rand('state',1);

  rate = 3/4; 
  framesize = 576;  

  mod_order = 4; 
  modulation = 'QPSK';
  mapping = 'gray';

  demod_type = 0;
  decoder_type = 0;
  max_iterations = 100;

  nbitspervocoderframe = 52;
  nvocoderframes = 8;
  nbitspermodemframe = 72;

  code_param = ldpc_init(rate, framesize, modulation, mod_order, mapping);

  data = [];
  r = []; 
  load interleaver.txt
  interleaver = interleaver + 1;

  % Encoder: Generate simulated vocoder data
  %          LPDC encode
  %          interleave           
  %          insert UW bits

  Nframes = 100;
  uw = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];

  % repeat same simulated vocoder data to ease testing

  vd = round( rand( 1, nbitspervocoderframe*nvocoderframes) );

  % pad data with zeros the size of UW

  vdpad = [vd zeros(1, length(uw))];

  % LDPC encode

  [codewordpad, s] = ldpc_enc(vdpad, code_param);
  code_param.code_bits_per_frame = length(codewordpad);
  code_param.symbols_per_frame = length(s);

  % remove padded zeros after encoding to leave room for UW bits (code
  % is systematic)

  codeword = [ codewordpad(1:length(vd)) codewordpad((length(vd)+length(uw)+1):length(codewordpad)) ];

  % interleave, insert UW bits, and pack bits (as C modulator likes packed bits)

  codeword_interleaved = interleave_bits(interleaver, codeword);
  codeword_interleaved_uw = [insert_uw(codeword_interleaved(1:length(vd)), uw) codeword_interleaved(length(vd)+1:length(codeword_interleaved)) ];
  packedcodeword = packmsb(codeword_interleaved_uw);

  cwfilename = strcat(filename,"_codeword.bin");
  fc=fopen(cwfilename,"wb");
  for nn = 1: Nframes        
      fwrite(fc,packedcodeword,"uchar");
  end
  fclose(fc);

  %printf("framesize: %d data_bits_per_frame: %d code_bits_per_frame: %d\n", ...
  %        framesize, code_param.data_bits_per_frame,  code_param.code_bits_per_frame);

  printf("Encoded %d LDPC codewords, saved in packed file: %s\n", Nframes, cwfilename);

  % Modulator: Modulate to QPSK symbols ------------------------------------------

  nbytespackedcodeword=length(packedcodeword);
  fc=fopen(cwfilename,"rb");
  mcwfilename = strcat(filename,"_modcodeword.bin");
  fm=fopen(mcwfilename,"wb");
  nbytespackedmodemframe = nbitspermodemframe/8;
  n = 0;

  [packedmodem, count] = fread(fc,nbytespackedmodemframe,"uchar");
  while (count == nbytespackedmodemframe)
      n++;
      unpackedmodem = unpackmsb(packedmodem);

      ii = 1;
      for i=1:2:length(unpackedmodem)
          mod_unpackedmodem(ii) = qpsk_mod(unpackedmodem(i:i+1));
          mod_unpackedmodem_float32(i) = real(mod_unpackedmodem(ii));
          mod_unpackedmodem_float32(i+1) = imag(mod_unpackedmodem(ii));
          ii += 1;
      end

      fwrite(fm, mod_unpackedmodem_float32, "float32");
      [packedmodem, count] = fread(fc,nbytespackedmodemframe,"uchar");
  end
  fclose(fc);
  fclose(fm);
  printf("Modulated %d modem frames to file: %s\n", n, mcwfilename);
endfunction


