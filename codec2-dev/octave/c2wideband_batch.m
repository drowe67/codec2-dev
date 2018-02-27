% c2wideband_batch.m
%
% Copyright David Rowe 2017
% This program is distributed under the terms of the GNU General Public License 
% Version 2

#{

  Octave script to batch process model parameters for wideband Codec 2.

  Outputs a set of model parameters that can be fed to c2sim for
  listening tests.  The companion c2wideband_fbf.m script is used to
  visualise the processing frame by frame
 
  c2sim -> dump files -> c2wideband_batch.m -> output model params -> c2sim -> play
 
  Usage:

    Build codec2 with -DDUMP (see codec2-dev/README), then generate dump files:

      ~/codec2-dev/build_linux/src$ ./c2sim ~/Desktop/c2_hd/speech_orig_16k.wav --Fs 16000 --dump speech

    Start Octave and generate the map file, this only needs to be done once:

      $ cd ~/codec2-dev/octave
      $ octave
      octave:1> c2wideband_batch("../build_linux/src/speech", "mode", "generate map")

    Then to run batch simulation and generate output speech:
   
       octave:1> c2wideband_batch("../build_linux/src/speech");

      ~/codec2-dev/build_linux/src$ ./c2sim ~/Desktop/c2_hd/speech_orig_16k.wav --Fs 16000 --phase0 --postfilter --amread speech_am.out --hmread speech_hm.out -o | play -t raw -r 16000 -s -2 -

    Testing C Port:

    $ cd codec2-dev/build_linux/unittest/tc2wideband
    $ ./tc2wideband
    $ cd ../../octave
    $ octave
    c2wideband_batch("../build_linux/src/speech", "verifyc", "../build_linux/unittest/tc2wideband_out.txt");

#}


function [surface mean_f] = c2wideband_batch(input_prefix, varargin)
  newamp;
  autotest;
  more off;

  max_amp = 160;
  mean_f = [];

  % defaults

  synth_phase = 1;
  output = 1;
  output_prefix = input_prefix;
  fit_order = 0;
  mode = "dct2";
  verifyc = 0;

  % load model paraneters
  
  model_name = strcat(input_prefix,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);

  % parse variable argument list

  if (length (varargin) > 0)

    ind = arg_exists(varargin, "mode");
    if ind
      mode =  varargin{ind+1};
    end

    % check for the "output_prefix" option

    ind = arg_exists(varargin, "output_prefix");
    if ind
      output_prefix =  varargin{ind+1};
    end

    ind = arg_exists(varargin, "no_output");
    if ind
      output = 0;
      synth_phase = 0;
    end

    ind = arg_exists(varargin, "verifyc");
    if ind
      verifyc = 1; output = 0; frames = 160;
      c_filename =  varargin{ind+1}
    end
  end

  if output && !strcmp(mode,"generate map")
    printf("output_prefix: %s\n",  output_prefix);
  end
  
  % Choose experiment to run test here -----------------------

  if strcmp(mode, "generate map")
    generate_map(model, K=30, "c2wideband_map");
    output = 0;
  end
  if strcmp(mode, "dct2")
    [model_ rate_K_surface] = experiment_rate_K_dct2(model(1:frames,:), 1);
    frames_out = rows(model_);
  end

  % optional verification against C port

  if verifyc
    load(c_filename);
    fg = 1;
    figure(fg++); clf; mesh(rate_K_surface); title('rate K surface');
    figure(fg++); clf; mesh(rate_K_surface_c); title('rate K surface C');
    figure(fg++); clf; mesh(rate_K_surface - rate_K_surface_c); title('difference');
  end
  
  % ----------------------------------------------------

  if output
    Am_out_name = sprintf("%s_am.out", output_prefix);
    fam  = fopen(Am_out_name,"wb"); 
 
    if synth_phase
      Hm_out_name = sprintf("%s_hm.out", output_prefix);
      fhm = fopen(Hm_out_name,"wb"); 
    end

    for f=1:frames_out
      %printf("%d ", f);   
      Wo = model_(f,1); L = min([model_(f,2) max_amp-1]); Am = model_(f,3:(L+2));
      %Wo = model(f,1); L = min([model(f,2) max_amp-1]); Am = model(f,3:(L+2));
      if Wo*L > pi
        printf("Problem: %d  Wo*L > pi\n", f);   
      end

      Am_ = zeros(1,max_amp); Am_(2:L) = Am(1:L-1); fwrite(fam, Am_, "float32");

      if synth_phase

        % synthesis phase spectra from magnitiude spectra using minimum phase techniques

        fft_enc = 256;
        phase = determine_phase(model_, f, fft_enc);
        assert(length(phase) == fft_enc);

        % sample phase at centre of each harmonic, not 1st entry Hm[1] in octave Hm[0] in C
        % is not used

        Hm = zeros(1, 2*max_amp);
        for m=1:L
          b = round(m*Wo*fft_enc/(2*pi));
          Hm(2*m) = cos(phase(b));
          Hm(2*m+1) = -sin(phase(b));
        end
        fwrite(fhm, Hm, "float32");    
      end
    end

    fclose(fam);
    if synth_phase
      fclose(fhm);
    end
  end % if output .....
 
  printf("\n")

endfunction
 

function ind = arg_exists(v, str) 
   ind = 0;
   for i=1:length(v)
      if strcmp(v{i}, str)
        ind = i;
      end     
    end
endfunction


% Create "map" from a training database.  The map tells us which order
% to read out and quantise DCT coeffs.  Could be approximated by a
% zig-zag pattern. 
%
% Example 3x3 map in a zig-zag pattern, so (1,1) is first coeff, (1,2)
% 2nd ....
% 
%  1 2 6 
%  2 5 7
%  4 8 9
%
% TODO: [ ] Come up with a better name than 
%       [ ] Script to convert this to a C header file.

function generate_map(model, K, map_filename)
  newamp; 
  c2wideband_const;

  [frames nc] = size(model);
  surface = resample_const_rate_f_mel(model, K, Fs);
  [nr nc] = size(surface);
  asurf = surface(1:dec:nr,:);
  [map rms_map mx mx_ind unwrapped_dcts] = create_map_rms(asurf, Nt, K);

  %printf("non zero coeffs: %d\n", sum(sum(map == 1)));
  figure(2); clf;
  mesh(map);
 
  printf("generating map file: %s\n", map_filename);
  save("-ascii", map_filename, "map");
endfunction


% ---------------------------------------------------------------------------------------
% rate K mel-resampling, high end correction, and DCT experiment workhorse

function [model_ rate_K_surface] = experiment_rate_K_dct2(model, plots=1)
  newamp;
  c2wideband_const;

  [frames nc] = size(model);

  % break into blocks of (Nt time samples) x (K freq samples)

  Nblocks = floor(frames/(Nt*dec));
  
  printf("frames: %d Nblocks: %d\n", frames, Nblocks);

  % map that defines order we read out and quantise DCT coeffs
  % TODO: for C port we need an Octave function to write Map to a C
  % include file

  map = load("c2wideband_map");

  % create arrays to reverse map quantiser_num to r,c Luts

  rmap = cmap = zeros(1,Nt*K);
  for r=1:Nt
    for c=1:K
       quantiser_num = map(r,c);
       rmap(quantiser_num) = r;
       cmap(quantiser_num) = c;
    end
  end

  % per-block processing ----------------------------------------------------

  % init a bunch of output variables

  rate_K_surface_ = zeros(Nblocks*Nt*dec, K);  
  sumnz = zeros(1,Nblocks);
  dct2_sd = zeros(1,Nblocks);
  model_ = [];

  for n=1:Nblocks
    st = (n-1)*dec*Nt+1; en = st + dec*Nt - 1;
    %printf("st: %d en: %d\n", st, en);

    [model_block_ adct2_sd qn rate_K_surface_block rate_K_surface_block_] = wideband_enc_dec(model(st:en,:), rmap, cmap);

    model_ = [model_; model_block_];

    % log these for plotting/development

    rate_K_surface(st:en,:) = rate_K_surface_block;
    rate_K_surface_(st:en,:) = rate_K_surface_block_;
    dct2_sd(n) = adct2_sd;
    sumnz(n) = qn;
  end

  fg = 1;
  
 figure(fg++); clf;  mesh(rate_K_surface); title('rate K surface block');
 
  if plots
    figure(fg++); clf; plot(sumnz); hold on; 
    plot([1 length(sumnz)],[mean(sumnz) mean(sumnz)]); hold off; title('Non Zero');
    figure(fg++); clf; plot(dct2_sd); title('DCT SD');
  end
  printf("average dct spectral distortion: %3.2f dB\n", mean(dct2_sd));
  printf("mean number of coeffs/DCT: %3.2f/%d\n", mean(sumnz), Nt*K);
  printf("coeffs/second: %3.2f\n", mean(sumnz)/(Nt*Tf*dec));
  printf("bits/s: %3.2f\n", 2.9*mean(sumnz)/(Nt*Tf*dec));

  % this measure just works on 20ms frames, not sure if that's correct

  dist = std((rate_K_surface_(1:dec:Nblocks*Nt*dec,:) - rate_K_surface(1:dec:Nblocks*Nt*dec,:))');
  
  if plots
    figure(fg++); clf; plot(dist); title('Rate K SD');
    printf("Rate K spectral distortion mean: %3.2f dB var: %3.2f\n", mean(dist), var(dist));
  end
endfunction


% Encode/decoder a 160ms block of model parameters
% TODO: (i) quantisation of DCT coeffs (ii) break into separate encoder and decoder functions

function [model_block_ dct2_sd qn rate_K_surface_block rate_K_surface_block_] = wideband_enc_dec(model_block, rmap, cmap)
    c2wideband_const;

    sim_quant = 1; % used to simulate quantisation, set to 1,2,4, etc
    dist_dB   = 2; % use enough coefficients to get this distortion ond DCT coeffs

    % Resample variable rate L vectors to fixed length rate K.  We have
    % left high end correction out for now, this is less of an issue
    % with a higher K

    [rate_K_surface_block rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model_block, K, Fs, 'para');

    % decimate down to 20ms time resolution, and DCT

    D = dct2(rate_K_surface_block(1:dec:Nt*dec,:));

    % So D is the 2D block of DCT coeffs at the encoder.  We want to
    % create a quantised version at the "decoder" E.  This loop copies
    % DCTs coeffs from D to E, until we get beneath a distortion
    % threshold.

    % This is essentially variable rate quantisation, but gives us
    % some idea of the final bit rate.  In practice we will also need
    % to constrain the total number of bits (ie bit rate), and
    % quantise each coefficient.

    % Turns out than mean SD (across many blocks/frames) is about the
    % same in the DCT domain as the rate K domain.  So we can just
    % measure MSE between D and E to estimate mean SD of the rate K
    % vectors after quantisation.

    E = mapped = zeros(Nt,K);

    qn = 0;
    adct2_sd = mean(std(D-E));
    while adct2_sd > dist_dB
      qn++;
      E(rmap(qn), cmap(qn)) = sim_quant*round(D(rmap(qn), cmap(qn))/sim_quant);
      adct2_sd = mean(std(D-E));
      %printf("qn %d %f\n", qn, adct2_sd);
    end

    % note neat trick to interpolate to 10ms frames despite dec to 20ms, this means
    % we don't need a separate decode side interpolator.

    dct2_sd = mean(std(D-E));
    rate_K_surface_block_ = idct2([sqrt(dec)*E; zeros(Nt*(dec-1), K)]);
    model_block_ = resample_rate_L(model_block, rate_K_surface_block_, rate_K_sample_freqs_kHz, Fs, 'para');
endfunction
