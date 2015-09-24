% ldpcdec.m
% David Rowe 31 Dec 2013
% 
% LDPC decoder test program, given a file of QPSK symbols (IQIQ floats), 
% performs frame sync, decodes, and measures BER.

function ldpcdec(filename, Eprob)

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
  if (nargin == 1)
    Eprob = 0.0;
  end

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

  Nc = 18;

  % repeat same simulated vocoder data to ease testing

  vd = round( rand( 1, nbitspervocoderframe*nvocoderframes) );

  % Decoder: Sync with LDPC frames, de-interleave, LDPC decode, strip off UW, measure BER -------

  mcwfilename = strcat(filename,"_modcodeword.bin");
  fm=fopen(mcwfilename,"rb");
  etfilename = strcat(filename,"_modcodeword_et.bin");
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

  Terrs = 0; Trawerrs = 0; Ferrs = 0; Tbits = 0; Tframes = 0; nerr = []; nrawerr = [];
  corr = []; n = 0;
  sync_state = 0; sync_count = 0;
  mod_unpackedmodem_log = [];
  sync_state_log = [];
  entered_track_log = [];
  sig_pwr_log = [];

  symbols = erasures = 0;

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
      %erasures = rand(1,length(mod_unpackedmodem)) < Eprob; 
      %mod_unpackedmodem(erasures) = 0;

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
 
          % determine BER stats of raw data before decoding

          raw_bits = zeros(1, code_param.data_bits_per_frame - length(uw));
          for i=1:(code_param.data_bits_per_frame - length(uw))/2
              raw_bits(2*(i-1)+1:2*i) = qpsk_demod(mod_codeword_deinter(i));
          end
          error_positions = xor(vd, raw_bits);
          Nerrs = sum(error_positions);
          Trawerrs += Nerrs;
          nrawerr(Tframes) = Nerrs;

          % Determine Es/N for each carrier. For this codeword we assume
          % across codeword (currently 320ms) signal is stationary.
          % So for each carrier signal level is constant, so we can
          % average across all symols of that carrier to get a better
          % estimate of the carrier power.  The spectral noise density
          % No will be the same for the bandwidth of each carrier.  So
          % we can use noise samples from all symbols together to get
          % a better estimate of the noise power.
          
          sig_pwr(Tframes,:) = zeros(1,Nc);
          noise_samples = [];
          for n=1:Nc

            % extract a vector of one carrier's symbols for this codeword
            % rotate so that decision boundaries are now real and imag axis

            r = mod_codeword(n:Nc:length(mod_codeword)) .* exp(j*pi/4);

            sig_est = mean(abs(r));

            % The noise is the variance of symbols (samples) about the actual symbol position
            % we reflect all symbols into the first quadrant to simplify things, as the actual
            % received symbol isn't matter, just the noise around it.  We model the received
            % symbol based on the estimated signal level.

            refl_symbols = abs(real(r)) + j*abs(imag(r));    
            est_symbols = exp(j*pi/4)*sig_est*ones(1,length(r));
            noise_samples = [ noise_samples (est_symbols - refl_symbols)];
                       
            sig_pwr(Tframes,n) = sig_est .^ 2;
          end
          noise_pwr(Tframes) = var(noise_samples);
          %plot(real(refl_symbols), imag(refl_symbols), '+');
          %hold on;
          %plot(real(exp(j*pi/4)*sig_est*ones(1,length(r))), imag(exp(j*pi/4)*sig_est*ones(1,length(r))), 'r+');
          %hold off;
          %printf("SNR: %f\n", 10*log10(sig_est*sig_est/noise_pwr(Tframes)));
 
          % Set erasures for carrier beneath a certain Es/N
          
          for n=1:Nc
              symbols++;
              EsN(n) = 10*log10(sig_pwr(Tframes,n)/noise_pwr(Tframes));
              if (EsN(n) < 1)
                 %mod_codeword(n:Nc:length(mod_codeword)) = 0;    
                 %printf("Tframes: %d n: %d EsN = %3.2fdB\n", Tframes, n, EsN(n)); 
                 erasures++;                                     
              end
          end

          % De-interleave again with erasures set ----------------------

          % remove UW symbols

          mod_codeword_no_uw = remove_uw_symbols(mod_codeword, code_param.data_bits_per_frame/2 - length(uw)/2, length(uw)/2);
          mod_codeword_no_uw = [mod_codeword_no_uw mod_codeword((code_param.data_bits_per_frame/2+1):code_param.code_bits_per_frame/2)];

          tmp = deinterleave_symbols(interleaver, mod_codeword_no_uw);
        
          % insert known symbols at end of data

          mod_codeword_deinter = [ tmp(1:(code_param.data_bits_per_frame/2 - length(uw)/2)) ...
                                 ones(1,length(uw)/2) * qpsk_mod([0 0]) ...
                                 tmp((code_param.data_bits_per_frame/2 - length(uw)/2+1):length(tmp)) ];

          % LDPC decode ------------------------------------------------

          detected_data = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, mod_codeword_deinter, EsNo);

          % unpack payload data, removing UW

          vd_rx = detected_data(1:(code_param.data_bits_per_frame - length(uw)));

          % measure coded BER

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

  printf("\nFrames: %d bits: %d errors: %d Raw BER = %f Coded BER = %f FER = %f\n", Tframes, Tbits, Terrs, Trawerrs/Tbits, Terrs/Tbits, Ferrs/Tframes);
  printf("Symbols: %d Erasures: %d  %f\n", symbols, erasures, erasures/symbols);
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
  plot(nrawerr);
  subplot(313)
  plot(nerr);

  figure(10);
  plot(10*log10(sig_pwr(:,3)./noise_pwr(:)),'b');
  hold on;
  plot(10+10*log10(noise_pwr(:)));
  plot(10+10*log10(sig_pwr(:,3)),'r');
%  for n=2:Nc
%    plot(n*10+10*log10(sig_pwr(:,n)./noise_pwr(:,n)));
%    plot(n*10+10*log10(sig_pwr(:,n)),'r');
%  end
  hold off;

  y = 1:Tframes;
  x = 1:Nc;
  z = 10*log10(sig_pwr(:,:)./((noise_pwr(:)*ones(1, Nc))));
  %printf("mean SNR = %3.2fdB\n", mean(z));
  figure(11);
  imagesc(x,y,z);
  figure(12);
  mesh(x,y,z);
  axis([1 Nc 1 Tframes 5 15]);

endfunction
