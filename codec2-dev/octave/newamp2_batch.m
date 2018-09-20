% newamp2_batch.m
%
% Copyright David Rowe 2018
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

#{

  Octave script to batch process model parameters using the new
  amplitude model, version 2.  Outputs another set of model parameters
  that can be fed to c2sim for listening tests.  The companion
  newamp2_fbf.m script is used to visualise the processing frame by frame
 
  c2sim -> dump files -> newamp1_batch.m -> output model params -> c2sim -> play
 
  Usage:

    build codec2 with -DDUMP - see codec2-dev/README, then:

    ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
    $ cd ~/codec2-dev/octave
    octave:14> newamp2_batch("../build_linux/src/hts1a")
    ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --pahw hts1a -o - | aplay -f S16

    1/ Bit stream development:

    $ ./c2sim ../../raw/vk5qi.raw --framelength_s 0.0125 --dump vk5qi_l --phase0 --lpc 10 --dump_pitch_e vk5qi_l_pitche.txt
    octave:526> newamp2_batch("../build_linux/src/vk5qi_l","output_prefix","../build_linux/src/vk5qi_l_dec", "mode", "encdec");
    $ ./c2sim ../../raw/vk5qi.raw --framelength_s 0.0125 --pahw vk5qi_l_dec --hand_voicing vk5qi_l_dec_v.txt -o - |  play -t raw -r 8000 -s -2 -

    2/ Generate a bit stream file:

    octave:101> newamp2_batch("../build_linux/src/vk5qi_l","no_output","mode", "enc", "bitstream", "vk5qi_2200_enc.c2");

    3/ Decode a bit stream file:

    octave:101> newamp2_batch("../build_linux/src/vk5qi_l","output_prefix","../build_linux/src/vk5qi_l_decbs", "mode", "dec", "bitstream", "vk5qi_2200_enc.c2");
    ./c2sim ../../raw/vk5qi.raw --framelength_s 0.0125 --pahw vk5qi_l_decbs --hand_voicing vk5qi_l_decbs_v.txt -o - |  play -t raw -r 8000 -s -2 -

   DCT based HQ/200 Candidate D (Sep 2018):

     1/ (h and s) trained:
     
       load train120_surf.mat
       newamp2; [s h] = design_huffman_enc(train120_surf, qstepdB=6, max_dcts=18);

     2/ Quantise all:

       all_surf_l = newamp2_batch("../build_linux/src/all_l", "output_prefix", "../build_linux/src/all_surf_l", "mode", "linear");
       newamp2; s_ = huffman_encode_surf(all_surf_l(:,2:35), qstepdB=6, max_dcts=20, max_bits=45, s, h);
       all_surf_l_ = all_surf_l; all_surf_l_(:,:) = 0; all_surf_l_(:,2:35)=s_;
        ewamp2_batch("../build_linux/src/all_l", "output_prefix", "../build_linux/src/all_surf_l_45a", "mode", "linear", "surf_in", all_surf_l_);
       $ ./c2sim ../../wav/all.wav --framelength_s 0.0125 --pahw all_surf_l_45a -o - | aplay -f S16_LE

     3/ Testing huffman encoder and decoder:
     
     newamp2; [s_ dc bits] = huffman_encode_surf(hts2a_surf_l(:,2:35), qstepdB=6, max_dcts=20, max_bits=45, s, h);
     newamp2; s_dec  = huffman_decode_surf(34, qstepdB=6, max_dcts=20, s, h, bits, dc);

#}


function surface = newamp2_batch(input_prefix, varargin)
  newamp2;
  more off;

  max_amp = 160; Fs = 8000;

  % defaults

  synth_phase = output = 1;
  output_prefix = input_prefix;
  mode = "linear";
  correct_rate_K_en = 0; quant = ""; vq_en = 0; M = 1; mask_en = 0;
  surface_in_en = 0;
  
  % parse variable argument list

  if (length (varargin) > 0)

    % check for the "output_prefix" option

    ind = arg_exists(varargin, "output_prefix");
    if ind
      output_prefix =  varargin{ind+1};
    end

    ind = arg_exists(varargin, "mode");
    if ind
      mode =  varargin{ind+1};
    end

    ind = arg_exists(varargin, "bitstream");
    if ind
      bitsream_filename =  varargin{ind+1};
    end

    ind = arg_exists(varargin, "no_output");
    if ind
      output = 0;
      synth_phase = 0;
    end

    if arg_exists(varargin, "deltaf")
      quant = "deltaf";
    end
    if arg_exists(varargin, "dct")
      quant = "dct";
    end
    if arg_exists(varargin, "dtlimit")
      quant = "dtlimit";
    end
    if arg_exists(varargin, "step")
      quant = "step";
    end
     
    ind = arg_exists(varargin, "vq");
    if ind
      quant = "vq";
      vq_en = 1; vq_name = varargin{ind+1}; vq_st = varargin{ind+2}; vq_en = varargin{ind+3};
    end

    ind = arg_exists(varargin, "M");
    if ind
      M = varargin{ind+1};
    end
    
    ind = arg_exists(varargin, "mask");
    if ind
      mask_en = 1;
      error_filename = varargin{ind+1};
    end

    correct_rate_K_en = arg_exists(varargin, "correct_rate_K");

    ind = arg_exists(varargin, "surf_in");
    if ind
      surface_in_en = 1;
      surface_in = varargin{ind+1};
    end
    
  end

  printf("output: %d", output);
  if (output)
    printf(" output_prefix: %s",  output_prefix);
  end
  printf(" mode: %s", mode);
  printf(" correct_rate_K: %d quant: %s vq_en: %d M: %d\n", correct_rate_K_en, quant, vq_en, M);
  if vq_en
    printf("vq_name: %s vq_st: %d vq_en: %d\n", vq_name, vq_st, vq_en);
  end
  if surface_in_en
    printf("surface_in_en: %d\n", surface_in_en);
  end
  
  model_name = strcat(input_prefix,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);
  printf("frames: %d\n", frames);
  
  voicing_name = strcat(input_prefix,"_pitche.txt");
  voicing = zeros(1,frames);
  
  if exist(voicing_name, "file") == 2
    printf("loading: %s\n", voicing_name);
    pitche = load(voicing_name);
    voicing = pitche(:, 3);
  end

  % Choose experiment to run test here -----------------------

  if (strcmp(mode,"linear") || strcmp(mode,"mel")) && !surface_in_en
    args.resampler = mode;

    vqs.en=vq_en; 
    if vq_en
      load(vq_name);
      vqs.st=vq_st;vqs.m=5;
      vqs.table = avq;
    end
    args.vq = vqs;
    args.quant = quant;
    args.correct_rate_K_en = correct_rate_K_en;
    args.M = M;
    args.plots = 1;
    [model_ surface sd_log] = experiment_resample(model, args);
  end

  # resample a linear surface, used to test deep learning ideas
  
  if strcmp(mode,"linear") && surface_in_en
    [model_ surface] = resample_surf(model, surface_in);
  end
  
  # combined bitstream encoded and decoder for development
  
  if strcmp(mode,"encdec")
    [bits rate_K_surface_] = candc_encoder(model, voicing);
    [model_ voicing_ rate_K_surface_dec_] = candc_decoder(bits);
    if exist("voicing_", "var")
      printf("voicing_ exists!\n");
    end
    % sanity check - should be a flat surface
    figure(7);
    en = 100;
    X = rate_K_surface_(1:en,:) - rate_K_surface_dec_(1:en,:);
    mesh(X(1:2:en,:))    
  end

  # stand alone model file -> bit stream file cand C encoder

  if strcmp(mode,"enc")
    [bits rate_K_surface_] = candc_encoder(model, voicing);
    fbit = fopen(bitsream_filename,"wb"); 
    fwrite(fbit, bits, "uchar");
    fclose(fbit);
  end

  # stand alone model file -> bit stream file cand D encoder

  if strcmp(mode,"encd")
    [bits rate_K_surface_] = candd_encoder(model, voicing);
    fbit = fopen(bitsream_filename,"wb"); 
    fwrite(fbit, bits, "uchar");
    fclose(fbit);
  end

  # stand alone bit stream file -> model file decoder

  if strcmp(mode,"dec")
    fbit = fopen(bitsream_filename,"rb"); 
    bits = fread(fbit, "uchar")';
    fclose(fbit);
    if mask_en
      % optional error masking, read in error file, count errors/frame
       errors_per_codec_frame = count_errors(error_filename);
    end
    [model_ voicing_ rate_K_surface_dec_] = candc_decoder(bits, errors_per_codec_frame);
    figure(1); mesh(rate_K_surface_dec_(1:240,:))
  end

  if strcmp(mode,"decd")
    fbit = fopen(bitsream_filename,"rb"); 
    bits = fread(fbit, "uchar")';
    fclose(fbit);
    errors_per_codec_frame = [];
    if mask_en
      % optional error masking, read in error file, count errors/frame
       errors_per_codec_frame = count_errors(error_filename);
    end
    [model_ voicing_ rate_K_surface_dec_] = candd_decoder(bits, errors_per_codec_frame);
    figure(1); mesh(rate_K_surface_dec_(1:240,:))
  end

  % ----------------------------------------------------

  if output
    Am_out_name = sprintf("%s_am.out", output_prefix);
    fam  = fopen(Am_out_name,"wb"); 
 
    Wo_out_name = sprintf("%s_Wo.out", output_prefix);
    fWo  = fopen(Wo_out_name,"wb"); 
  
    if synth_phase
      Hm_out_name = sprintf("%s_hm.out", output_prefix);
      fhm = fopen(Hm_out_name,"wb"); 
    end

    for f=1:frames
      printf("%d\r", f);   
      Wo = model_(f,1); L = min([model_(f,2) max_amp-1]); Am = model_(f,3:(L+3));
      if Wo*L > pi
        printf("Problem: %d  Wo*L > pi Wo: %f F0: %f L: %d Wo*L: %f\n", f, Wo, Wo*Fs/(2*pi), L, Wo*L);   
      end

      Am_ = zeros(1,max_amp); Am_(2:L+1) = Am(1:L); fwrite(fam, Am_, "float32");
      fwrite(fWo, Wo, "float32");

      if synth_phase

        % synthesis phase spectra from magnitiude spectra using minimum phase techniques

        fft_enc = 128;
        phase = determine_phase(model_, f, fft_enc);
        assert(length(phase) == fft_enc);

        % sample phase at centre of each harmonic, not 1st entry Hm(1:2) in octave Hm[0] in C
        % is not used

        Hm = zeros(1, 2*max_amp);
        for m=1:L
          b = round(m*Wo*fft_enc/(2*pi));
          Hm(2*m+1) = cos(phase(b));
          Hm(2*m+2) = sin(phase(b));
        end
        fwrite(fhm, Hm, "float32");
      end
    end
 
   fclose(fam);
   fclose(fWo);
   if synth_phase
     fclose(fhm);
   end

   % save voicing file
  
   if exist("voicing_", "var") && output
     v_out_name = sprintf("%s_v.txt", output_prefix);
     printf("writing: %s\n", v_out_name);
     fv  = fopen(v_out_name,"wt"); 
     for f=1:length(voicing_)
       fprintf(fv,"%d\n", voicing_(f));
     end
     fclose(fv);
    end
  end

  printf("\n")

endfunction
 

function ind = arg_exists(v, str) 
   ind = 0;
   for i=1:length(v)
      if !ind && strcmp(v{i}, str)
        ind = i;
      end     
    end
endfunction


% Basic unquantised rate K sampling (mel or linear) then back to rate L

function [model_ rate_K_surface_ sd_log delta_K] = experiment_resample(model, args)
  [frames nc] = size(model);
  K = 20; Fs = 8000;
  AmdB = zeros(frames, 160);
  melvq;

  resampler = args.resampler;
  quant = args.quant;
  vq    = args.vq;
  plots = args.plots;
  correct_rate_K_en = args.correct_rate_K_en;
  M = args.M;
  
  if strcmp(resampler, "mel")
    K = 20;
  end
  if strcmp(resampler, "linear")
    K = 40; step = (Fs/2000)/K;
    rate_K_sample_freqs_kHz = [0.1:step:4];
    K = length(rate_K_sample_freqs_kHz);
  end
  
  rate_K_surface = rate_K_surface_ = zeros(frames, K);
    
  for f=1:frames
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    energy_L(f) = sum(Am.^2);
    AmdB(f,1:L) = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*Fs/(2000*pi);
    if strcmp(resampler, "mel")
      [rate_K_vec rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model(f,:), K, Fs, 'lanc');
    end
    if strcmp(resampler, "linear")
      rate_K_vec = resample_const_rate_f(model(f,:), rate_K_sample_freqs_kHz, Fs, 'lanc');
    end
    rate_K_vec_lin = 10 .^ (rate_K_vec/20);
    energy_K(f) = sum(rate_K_vec_lin .^ 2);
    #g =  sqrt(energy_L(f)./energy_K(f));
    #rate_K_vec_lin *= g;
    energy_K(f) = sum(rate_K_vec_lin .^ 2);
    
    if correct_rate_K_en
      [tmp_ AmdB_] = resample_rate_L(model(f,:), rate_K_vec, rate_K_sample_freqs_kHz, Fs, 'lancmel');
      [rate_K_vec_corrected orig_error error nasty_error_log nasty_error_m_log] = correct_rate_K_vec(rate_K_vec, rate_K_sample_freqs_kHz, AmdB, AmdB_, K, Wo, L, Fs);
      rate_K_surface(f,:) = rate_K_vec_corrected;
    else
      rate_K_surface(f,:) = rate_K_vec;
    end
  end

  % Whole thing is quantised to 6dB steps, as that doesn't seem to
  % introduce much distortion
  
  if strcmp(quant,"step");
    rate_K_surface_ = 6*round(rate_K_surface/6);
  end
  
  % experiment to limit frame by frame changes, this will make quantisation easier

  if strcmp(quant,"dtlimit");
    rate_K_surface_prev_ = zeros(1,K);
    delta_K_log = zeros(frames,K);

    for f=1:M:frames

      % find delta

      delta_K = rate_K_surface(f,:) - rate_K_surface_prev_;

      % limit to -6,0,+6 dB 

      delta_K = 6*round(delta_K/6);
      delta_K = min(delta_K, +12);
      delta_K = max(delta_K, -12);

      % optional VQ

      if nargin == 5
        [res output_vec ind] = mbest(vq.table, delta_K(:,vq.st:vq.en), vq.m);
        delta_K(:,vq.st:vq.en) = output_vec;
      end
      
      rate_K_surface_prev_ += delta_K;
      rate_K_surface_(f,:) = rate_K_surface_prev_;
      delta_K_log(f,:) = delta_K;
    end
  end

  % encoding of delata_f, also limits delta_f changes, (hopefully) making quantisation easier

  E_log = [];
  if strcmp(quant,"deltaf");
    for f=1:M:frames

      % used to gather stats for frame energy
      
      E = 6*round(rate_K_surface(f,3)/6);
      E_log = [E_log E/6];

      % Determine frame energy as k=3 sample and quantise From looking
      % at values from all.wav it looks like 4 bits is enough

      levels = (0:15)*6;
      [E_ E_index] = quantise(levels, rate_K_surface(f,3));

      [rate_K_surface_(f,:) spec_mag_bits] = deltaf_quantise_rate_K_huff(rate_K_surface(f,:), E_, 44);

      # test of decoder
      
      dec_rate_K_surface_ = deltaf_decode_rate_K_huff(spec_mag_bits, E_, K, 44);
      assert (rate_K_surface_(f,:) == dec_rate_K_surface_);

      # quantise pitch

      Wo = model(f,1);
      Wo_index = encode_log_Wo(Wo, 7);

      #{
        TODO:
             [ ] save to bits file
             [ ] reconsruct frame here from indexes
             [ ] compare results
             [ ] then refactor into sep enc and dec functions,make sure identical results
      #}
      
      % printf("length: %d\n", length(bits));

      #{
      % add some experimental random noise of +/- one entry to simulate discrete gain-shape
      % VQ errors.  Make sure we choose a different position for each error.  used during VQ
      % development that didn't quite make it ....

      kmax=20;
      krecord = zeros(1,kmax);
      for i=1:4
        k = ceil(kmax*rand(1,1));
        while krecord(k)
          k = ceil(kmax*rand(1,1));
        end
        val = 6-12*floor(2*rand(1,1));
        krecord(k) = val;
        rate_K_surface_(f,k) += val;
        %printf("f: %d k: %d val: %d\n", f, k, val);
      end
      krecord
      #}
    end
  end

  if strcmp(quant,"dct");
    for f=1:M:frames
      rate_K_surface_(f,:) = dct_quantise_rate_K(rate_K_surface(f,:));
    end
  end 
  if strcmp(quant,"");
    rate_K_surface_ = rate_K_surface;
  end
  
  if strcmp(quant,"vq")
    % use orig for samples we don't quantise
    rate_K_surface_ = rate_K_surface;

    size(vq.table)
    target = rate_K_surface(:,vq.st:vq.en);
    [res output_vecs ind] = mbest(vq.table, target, vq.m);
    rate_K_surface_(:,vq.st:vq.en) = output_vecs;
  end

  % optional interpolation, if decimation M > 1

  if M > 1
    for f=1:M:frames-M
      left_f = f; right_f = f+M;
      %printf("%d %d\n", left_f, right_f);
      left_vec = rate_K_surface_(left_f, :); right_vec = rate_K_surface_(right_f, :);
      sample_points = [left_f right_f];
      resample_points = left_f+1:right_f-1;
      for k=1:K
        rate_K_surface_(resample_points, k) = interp_linear(sample_points, [left_vec(k) right_vec(k)], resample_points);
      end
    end
    rate_K_surface_(frames-M+1:frames,:) = rate_K_surface_(frames-M,1);    
  end

  if plots
    en = min(100, length(rate_K_surface_));
    st = en-49;
    figure(1); clf; mesh(rate_K_surface_(st:en,:)); ylabel('time'); xlabel('Mel freq');
  end

  % back to rate L
  
  if strcmp(resampler, "mel")
    [model_ AmdB_ ] = resample_rate_L(model, rate_K_surface_, rate_K_sample_freqs_kHz, Fs, 'lancmel');
  end
  if strcmp(resampler, "linear")
    [model_ AmdB_ ] = resample_rate_L(model, rate_K_surface_, rate_K_sample_freqs_kHz, Fs, 'lanc');
  end
  for f=1:frames
    Wo = model_(f,1);
    L = model_(f,2);
    Am_ = model_(f,3:(L+2));
    energy_L_(f) = sum(Am_.^2);
  end
  
  % calculate SD

  sd_log = []; energy_log = [];
  for f=1:frames
    L = model(f,2);
    asd = std(AmdB(f,1:L) - AmdB_(f,1:L));
    energy_log = [energy_log mean(AmdB(f,1:L))];
    sd_log = [sd_log asd];
  end

  if plots
    figure(2); clf; l = length(sd_log); ax=plotyy((1:l),sd_log,(1:l),energy_log);
    title('SD againstframe'); xlabel('Frame'); ylabel (ax(1), "SD"); ylabel (ax(2), "Energy");

    figure(3); clf; plot(hist(sd_log)); title('SD histogram');
    figure(4); clf; plot(energy_log, sd_log, '+'); title('Scatter of SD against energy');
    xlabel('Energy'); ylabel('SD'); axis([-10 60 0 7]);

    figure(5);
    plot(10*log10(energy_L)); hold on;
    plot(10*log10(energy_K),"g"); plot(10*log10(energy_L_),"r+"); hold off;
  end
  
  % return delta_K for training
  
  if strcmp(quant,"dtlimit")
    rate_K_surface_ = delta_K_log;
    l = min(100, length(rate_K_surface_));
    figure(5); clf; mesh(rate_K_surface_(1:l,:));
  end

  % this code useful to explore outliers, adjust threshold based on plot(sd_log)

  plot_outliers = 0;
  if plot_outliers && (asd > 7)
    for f=1:frames
      L = model(f,2);
      figure; plot(AmdB(f,1:L), 'b+-'); hold on; plot(AmdB_(f,1:L), 'r+-'); hold off;
    end
  end
  printf("mean SD %4.2f\n", mean(sd_log));
  printf("Emin: %f Emax: %f\n", min(E_log), max(E_log));

  #{
  % Plot E stats, by viewing these decided on 4 bit/frame for energy

  if length(E_log)
    figure(5);
    plot(E_log)
    figure(6);
    hist(E_log);
  end
  #}
endfunction


% Takes surface as input, back to model_, for deep learning experiments

function [model_ rate_K_surface_ sd_log delta_K] = resample_surf(model, rate_K_surface_)
  Fs = 8000;

  K = 40; step = (Fs/2000)/K;
  rate_K_sample_freqs_kHz = [0.1:step:4];
  K = length(rate_K_sample_freqs_kHz);
  
  % back to rate L
  
  [model_ AmdB_ ] = resample_rate_L(model, rate_K_surface_, rate_K_sample_freqs_kHz, Fs, 'lanc');
endfunction


% Candidate C, model to bit-stream encoder

function [bits rate_K_surface_] = candc_encoder(model, voicing)
  newamp2_const;
  [frames nc] = size(model);
  AmdB = zeros(frames, max_amp);
  bits = [];
  
  rate_K_surface = rate_K_surface_ = zeros(frames, K);
    
  for f=1:M:frames
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    AmdB(f,1:L) = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*Fs/(2000*pi);
    [rate_K_vec rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model(f,:), K, Fs, 'lanc');
    rate_K_surface(f,:) = rate_K_vec;
  end

  for f=1:M:frames

    % Determine frame energy as k=3 sample and quantise From looking
    % at values from all.wav it looks like 4 bits is enough

    levels = (0:15)*6;
    [E_ E_index] = quantise(levels, rate_K_surface(f,3));
    E_bits = index_to_bits(E_index-1, 4);
    [rate_K_surface_(f,:) spec_mag_bits] = deltaf_quantise_rate_K_fixed(rate_K_surface(f,:), E_, 44);

    # quantise pitch

    Wo = model(f,1);
    Wo_index = encode_log_Wo(Wo, 7);
    Wo_bits = index_to_bits(Wo_index, 7);

    # pack bits

    bits_frame = [Wo_bits E_bits voicing(f) spec_mag_bits ];
    bits = [bits bits_frame];

    if f < 10
      printf("f: %d Wo: %2d E: %d %2.0f v: %d\n", f, Wo_index, E_index, E_, voicing(f));
    end
  end

endfunction


% Candidate C, bit stream to model decoder

function [model_ voicing rate_K_surface_] = candc_decoder(bits, errors_per_codec_frame=[], ber = 0.0)
  newamp2_const;
  rows = floor(length(bits)/Nbitspercodecframe);
  frames = rows*M;
  
  rate_K_surface_ = zeros(frames, K);
  model_ = zeros(frames, max_amp+2);
  voicing = zeros(1,frames);
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K, 100, 0.95*Fs/2);
  Tbits = Terrs = 0;
  abits = zeros(1, Nbitspercodecframe);
  Nerrs = 0; av_level = 1;

  level_log = level_adj_log = [];
  
  error_thresh = 0;
  
  r = 1;
  for f=1:M:frames
    abits = bits(1:Nbitspercodecframe);
    bits = bits(Nbitspercodecframe+1:end);

    % optional insertion of bit errors for testing
    
    if ber > 0.0
      [abits(13:56) nerr] = insert_bit_error(abits(13:56), 0.00);
      Terrs += nerr;
      Tbits += 45;
    end
    
    % extract information from bit stream

    spec_mag_bits = abits(13:Nbitspercodecframe);
    E_bits = abits(8:11);
    Wo_bits = abits((1:7));
    voicing(f) = abits(12);
    Wo_index = bits_to_index(Wo_bits, 7);
    Wo_ = decode_log_Wo(Wo_index, 7); L_ = floor(pi/Wo_);
    E_index = bits_to_index(E_bits, 4);
    E_ = E_index*6;

    #{
    if f < 10
      printf("f: %d Wo: %2d E: %d %2.0f v: %d\n", f, Wo_index, E_index+1, E_, voicing(f));
    end
    #}
    
    % decode into rate K vec
    
    rate_K_surface_(f,:) = deltaf_decode_rate_K_fixed(spec_mag_bits, E_, K, 44);

    level = max(rate_K_surface_(f,:));
    level_log = [level_log level];
    
    if length(errors_per_codec_frame) >= r
      Nerrs = errors_per_codec_frame(r);
    end

    if Nerrs > error_thresh
      % if errors, want to avoid loud cracks, but we also don't want to simply mute.  So
      % adjust level to match recent average, and let that decay if long stream of error
      % frames to gradually mute.

      adjustment = av_level/level;
      if adjustment < 1
        rate_K_surface_(f,:) *= adjustment;
        printf("f: %3d r: %3d e: %2d level: %3.1f av_level: %3.1f adjust: %3.2f\n", f,r,Nerrs, level, av_level, adjustment);
      end
      av_level = av_level*0.9;      
    else
      % update average level
      av_level = av_level*0.9 + level*0.1;
      printf("f: %3d r: %3d e: %2d level: %3.1f av_level: %3.1f\n", f,r,Nerrs, level, av_level);
    end
       
    level = max(rate_K_surface_(f,:));
    level_adj_log = [level_adj_log level];
    
    model_(f,1) = Wo_; model_(f,2) = L_;
    r++;
  end

  if length(errors_per_codec_frame)
    figure(3); clf; nplot = 80*3;
    subplot(211,"position",[0.1 0.8 0.8 0.15]);
    plot(errors_per_codec_frame(1:nplot))
    subplot(212,"position",[0.1 0.05 0.8 0.7]);
    plot(level_log(1:nplot),'b'); hold on; plot(level_adj_log(1:nplot),'g+'); hold off;
  end
  
  if ber > 0.0
    printf("Tbits: %d Terrs: %d BER: %4.3f\n", Tbits, Terrs, Terrs/Tbits);
  end
  
  % interpolation in time
    
  for f=1:M:frames-M

    % interpolate rate K ampl samples
    
    left_f = f; right_f = f+M;
    left_vec = rate_K_surface_(left_f, :); right_vec = rate_K_surface_(right_f, :);
    sample_points = [left_f right_f];
    assert(M==2);
    resample_point = left_f+1;
    for k=1:K
      rate_K_surface_(resample_point, k) = interp_linear(sample_points, [left_vec(k) right_vec(k)], resample_point);
    end

    % interpolate Wo and voicing

    Wo_ = 2*pi/160; v = 0;
    if voicing(left_f) && voicing(right_f)
      Wo_ = (model_(left_f,1) + model_(right_f,1))/2;
      v = 1;
    end
    L_ = floor(pi/Wo_);
    model_(resample_point,1) = Wo_; model_(resample_point,2) = L_;
    voicing(resample_point) = v;
  end
  rate_K_surface_(frames-M+1:frames,:) = rate_K_surface_(frames-M,1);    
  model_(frames-M+1:frames,1) = model_(frames-M,1);    
  model_(frames-M+1:frames,2) = model_(frames-M,2);    

  % back to rate L amplitude samples
  
  model_ = resample_rate_L(model_, rate_K_surface_, rate_K_sample_freqs_kHz, Fs, 'lancmel');

endfunction


function [bits nerrs] = insert_bit_error(bits, ber)
  newamp2_const;
  p = rand(1,length(bits));
  error_mask = p < ber;
  bits = xor(bits, error_mask);
  nerrs = sum(error_mask);
endfunction


% Counts errors in protected bits.  Simulates failure of FEC to decode, which
% we can detect

function errors_per_codec_frame = count_errors(error_filename)
  newamp2_const;
  ferr = fopen(error_filename,"rb"); 
  errors = fread(ferr, "uchar")';
  fclose(ferr);
  frames = floor(length(errors)/Nbitspercodecframe)
  errors_per_codec_frame = zeros(1,frames);
  for f=1:frames
    st = (f-1)*Nbitspercodecframe + 1; en = st + Nprotectedbitspercodecframe - 1;
    errors_per_codec_frame(f) = sum(errors(st:en));
  end
  figure(2);
  plot(errors_per_codec_frame);
endfunction


% Candidate D (Huffman encoded DCTs), model to bit-stream encoder.  Note rate_K_surface
% is returned without energy being quantised, so not quite the same as decoder output.

function [bits rate_K_surface_] = candd_encoder(model, voicing)
  newamp2_candd_const;
  [frames nc] = size(model);
  AmdB = zeros(frames, max_amp);
  bits = [];
  load huffman.mat;
  
  rate_K_surface = rate_K_surface_ = zeros(frames, K);

  % extract vectors from model file and convert rate L amplitude samples to fixed rate K
  
  step = (Fs/2000)/K;
  rate_K_sample_freqs_kHz = [0.1:step:4];

  for f=1:frames
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    AmdB(f,1:L) = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*Fs/(2000*pi);
    [rate_K_vec rate_K_sample_freqs_kHz] = resample_const_rate_f(model(f,:), rate_K_sample_freqs_kHz, Fs, 'lanc');
    rate_K_surface(f,:) = rate_K_vec;
  end

  % huffman encode rate K amplitude samples

  [s_ dc spec_mag_bits] = huffman_encode_surf(rate_K_surface(:,2:35), qstepdB=6, max_dcts=20, max_bits=44, s, h);
  rate_K_surface_(:, 2:35) = s_;
  
  for f=1:dec:frames

    % Determine frame energy from DCT DC value. Looking at values from
    % all.wav it looks like 4 bits over 0 to 60dB is enough

    levels = (0:15)*4;
    [E_ E_index] = quantise(levels, dc(f));
    E_bits = index_to_bits(E_index-1, 4);

    # quantise pitch

    Wo = model(f,1);
    Wo_index = encode_log_Wo(Wo, 7);
    Wo_bits = index_to_bits(Wo_index, 7);

    # pack bits

    bits_frame = [Wo_bits E_bits voicing(f) spec_mag_bits{f}];
    bits = [bits bits_frame];

    if (f > 20) && (f < 30)
      printf("f: %3d Wo: %2d E: %2d %2.0f v: %d\n", f, Wo_index, E_index, E_, voicing(f));
    end
  end

  figure(1); clf; mesh(rate_K_surface_);
endfunction


% Candidate D, bit stream to model decoder

function [model_ voicing rate_K_surface_] = candd_decoder(bits, errors_per_codec_frame=[], ber = 0.0)
  newamp2_candd_const;
  load huffman.mat;
  rows = floor(length(bits)/Nbitspercodecframe);
  frames = rows*dec;
  
  rate_K_surface_ = zeros(frames, K);
  model_ = zeros(frames, max_amp+2);
  voicing = zeros(1,frames);
  step = (Fs/2000)/K;
  rate_K_sample_freqs_kHz = [0.1:step:4];
  Tbits = Terrs = 0;
  abits = zeros(1, Nbitspercodecframe);
  Nerrs = 0; av_level = 1;
  spec_mag_bits  = cell(frames,1);
  level_log = level_adj_log = [];
  E = zeros(1,frames);
  
  error_thresh = 0;
  
  r = 1;
  for f=1:dec:frames
    abits = bits(1:Nbitspercodecframe);
    bits = bits(Nbitspercodecframe+1:end);

    % optional insertion of bit errors for testing
    
    if ber > 0.0
      [abits(13:56) nerr] = insert_bit_error(abits(13:56), 0.00);
      Terrs += nerr;
      Tbits += 45;
    end
    
    % extract information from bit stream

    spec_mag_bits{f} = abits(13:Nbitspercodecframe);
    E_bits = abits(8:11);
    Wo_bits = abits((1:7));
    voicing(f) = abits(12);
    Wo_index = bits_to_index(Wo_bits, 7);
    Wo_ = decode_log_Wo(Wo_index, 7); L_ = floor(pi/Wo_);
    E_index = bits_to_index(E_bits, 4);
    E_(f) = E_index*4;

    if (f > 20) && (f < 30)
      printf("f: %d Wo: %2d E: %d %2.0f v: %d\n", f, Wo_index, E_index+1, E_(f), voicing(f));
    end
    
    model_(f,1) = Wo_; model_(f,2) = L_;
    r++;
  end

  % decode into rate K vec
  % TODO: make this frame by frame and run in loop above

  s_ = huffman_decode_surf(K=34, qstepdB=6, max_dcts=20, s, h, spec_mag_bits, E_);
  rate_K_surface_(:,2:35) = s_;

  if length(errors_per_codec_frame)
    figure(3); clf; nplot = 80*3;
    subplot(211,"position",[0.1 0.8 0.8 0.15]);
    plot(errors_per_codec_frame(1:nplot))
    subplot(212,"position",[0.1 0.05 0.8 0.7]);
    plot(level_log(1:nplot),'b'); hold on; plot(level_adj_log(1:nplot),'g+'); hold off;
  end
  
  if ber > 0.0
    printf("Tbits: %d Terrs: %d BER: %4.3f\n", Tbits, Terrs, Terrs/Tbits);
  end
  
  % interpolation in time
    
  for f=1:dec:frames-dec

    % interpolate rate K ampl samples
    
    left_f = f; right_f = f+dec;
    sample_points = [left_f right_f];
    assert(dec==2);
    resample_point = left_f+1;

    % interpolate Wo and voicing

    Wo_ = 2*pi/160; v = 0;
    if voicing(left_f) && voicing(right_f)
      Wo_ = (model_(left_f,1) + model_(right_f,1))/2;
      v = 1;
    end
    L_ = floor(pi/Wo_);
    model_(resample_point,1) = Wo_; model_(resample_point,2) = L_;
    voicing(resample_point) = v;
  end
  rate_K_surface_(frames-dec+1:frames,:) = rate_K_surface_(frames-dec,1);    
  model_(frames-dec+1:frames,1) = model_(frames-dec,1);    
  model_(frames-dec+1:frames,2) = model_(frames-dec,2);    

  % back to rate L amplitude samples
  
  model_ = resample_rate_L(model_, rate_K_surface_, rate_K_sample_freqs_kHz, Fs, 'lanc');

  figure(2); clf; mesh(rate_K_surface_);
endfunction
