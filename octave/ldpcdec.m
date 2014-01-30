% ldpcdec.m
% David Rowe 31 Dec 2013
% 
% LDPC decoder test program, given a file of QPSK symbols (IQIQ floats), 
% performs frame sync, decodes, and measures BER.

function ldpcdec(filename)

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
  EsNo = 4;
  Eprob = 0.0;

  nbitspervocoderframe = 52;
  nvocoderframes = 8;
  nbitspermodemframe = 72;

  code_param = ldpc_init(rate, framesize, modulation, mod_order, mapping);
  code_param.code_bits_per_frame = 576;

  data = [];
  r = []; 
  load interleaver.txt
  interleaver = interleaver + 1;

  Nframes = 100;
  uw = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];

  % repeat same simulated vocoder data to ease testing

  vd = round( rand( 1, nbitspervocoderframe*nvocoderframes) );

  % Decoder: Sync with LDPC frames, de-interleave, LDPC decode, strip off UW, measure BER -------

  mcwfilename = strcat(filename,"_modcodeword.bin");
  fm=fopen(mcwfilename,"rb");
  etfilename = strcat(filename,"_et.bin");
  fet=fopen(etfilename ,"rb");
  epfilename = strcat(filename,".err");
  fep=fopen(epfilename ,"wb");
  printf("Input QPSK symbols: %s\n", mcwfilename);
  if (fet == -1)
    printf("Input entered track file: none\n");
  else
    printf("Input entered track file: %s\n", etfilename);
  end
  printf("Output error pattern file: %s\n", epfilename);

  mod_uw = build_mod_uw(uw, 2*length(vd)/length(uw));

  mod_codeword = zeros(1, code_param.code_bits_per_frame/2);
  lmod_codeword = code_param.code_bits_per_frame/2;

  Terrs = 0; Ferrs = 0; Tbits = 0; Tframes = 0; nerr = [];
  corr = []; n = 0;
  sync_state = 0; sync_count = 0;
  mod_unpackedmodem_log = [];
  sync_state_log = [];
  entered_track_log = [];

  [mod_unpackedmodem_float32, count] = fread(fm,nbitspermodemframe, "float32");
  if (fet == -1)
      entered_track = 0;
  else
     entered_track = fread(fet, 1, "int");
  end

  while (count == nbitspermodemframe)
      n++;

      mod_unpackedmodem = mod_unpackedmodem_float32(1:2:nbitspermodemframe) + j*mod_unpackedmodem_float32(2:2:nbitspermodemframe);
      mod_unpackedmodem_log = [mod_unpackedmodem_log mod_unpackedmodem];
      erasures = rand(1,length(mod_unpackedmodem)) < Eprob; 
      mod_unpackedmodem(erasures) = 0;

      % keep buffer of one entire codeword

      mod_codeword(1:lmod_codeword-length(mod_unpackedmodem)) = mod_codeword(length(mod_unpackedmodem)+1:lmod_codeword);
      mod_codeword(lmod_codeword-length(mod_unpackedmodem)+1:lmod_codeword) = mod_unpackedmodem;

      [uw_sync corr(n)] = look_for_uw(mod_codeword(1:length(mod_uw)), mod_uw);

      next_sync_state = sync_state;
      if ((sync_state == 0) && (uw_sync == 1))
        next_sync_state = 1;
        sync_count = 0;
      end
      if ((sync_state == 1) && (entered_track != 0))
        next_sync_state = 0;
      end
      sync_state = next_sync_state;
      sync_state_log = [sync_state_log sync_state];
      entered_track_log = [entered_track_log entered_track];

      if (sync_state && (sync_count == 0))
          Tframes++;

          % remove UW symbols

          mod_codeword_no_uw = remove_uw_symbols(mod_codeword, code_param.data_bits_per_frame/2 - length(uw)/2, length(uw)/2);
          mod_codeword_no_uw = [mod_codeword_no_uw mod_codeword((code_param.data_bits_per_frame/2+1):code_param.code_bits_per_frame/2)];

          % de-interleave

          tmp = deinterleave_symbols(interleaver, mod_codeword_no_uw);
        
          % insert known symbols at end of data

          mod_codeword_deinter = [ tmp(1:(code_param.data_bits_per_frame/2 - length(uw)/2)) ...
                                 ones(1,length(uw)/2) * qpsk_mod([0 0]) ...
                                 tmp((code_param.data_bits_per_frame/2 - length(uw)/2+1):length(tmp)) ];
 
          % LDPC decode

          detected_data = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, mod_codeword_deinter, EsNo);

          % unpack payload data, removing UW

          vd_rx = detected_data(1:(code_param.data_bits_per_frame - length(uw)));

          % measure BER

          error_positions = xor(vd, vd_rx);
          Nerrs = sum(error_positions);
          if Nerrs>0, fprintf(1,'x'); Ferrs++; ,  else fprintf(1,'.'),  end
          Tbits += length(vd);
          Terrs += Nerrs;
          nerr(Tframes) = Nerrs;

          % save error patterns is simulated vocoder data to disk

          fwrite(fep, error_positions, "short");
          
      end

      if (sync_state)
          sync_count++;
          if (sync_count == 8)
              sync_count = 0;
          end
      end

      % read in one modulated modem frame at a time

      [mod_unpackedmodem_float32, count] = fread(fm, nbitspermodemframe, "float32");
      if (fet == -1)
          entered_track = 0;
      else
          entered_track = fread(fet, 1, "int");
      end
  end

  fclose(fep);

  fprintf(1,"\nFrames: %d bits: %d errors: %d BER = %f FER = %f\n", Tframes, Tbits, Terrs, Terrs/Tbits, Ferrs/Tframes);

  figure(8)
  clf;
  [n m] = size(mod_unpackedmodem_log);
  plot( real(mod_unpackedmodem_log), imag(mod_unpackedmodem_log), '+')
  axis([-2 2 -2 2]);
  title('Scatter Diagram');

  figure(9)
  subplot(311)
  plot(sync_state_log);
  subplot(312)
  plot(entered_track_log);
  subplot(313)
  plot(nerr);
endfunction
