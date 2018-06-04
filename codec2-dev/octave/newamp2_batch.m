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
    ~/codec2-dev/build_linux/src$ /c2sim ../../raw/hts1a.raw --pahw hts1a -o - | aplay -f S16
 
#}


function surface = newamp2_batch(input_prefix, varargin)
  newamp2;
  more off;

  max_amp = 160;

  % defaults

  synth_phase = output = 1;
  output_prefix = input_prefix;
  mode = "mel";
  
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

    ind = arg_exists(varargin, "no_output");
    if ind
      output = 0;
      synth_phase = 0;
    end
  end

  printf("output: %d\n", output);
  if (output)
    printf("output_prefix: %s\n",  output_prefix);
  end
  printf("mode: %s\n", mode);
  
  model_name = strcat(input_prefix,"_model.txt");
  model = load(model_name);
  [frames nc] = size(model);

  voicing_name = strcat(input_prefix,"_pitche.txt");
  voicing = zeros(1,frames);
  
  if exist(voicing_name, "file") == 2
    pitche = load(voicing_name);
    voicing = pitche(:, 3);
  end

  % Choose experiment to run test here -----------------------

  if strcmp(mode, 'mel')
    [model_ surface sd_log] = experiment_mel_freq(model, 0, 1, voicing);    
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
      %printf("%d ", f);   
      Wo = model_(f,1); L = min([model_(f,2) max_amp-1]); Am = model_(f,3:(L+2));
      if Wo*L > pi
        printf("Problem: %d  Wo*L > pi\n", f);   
      end

      Am_ = zeros(1,max_amp); Am_(2:L) = Am(1:L-1); fwrite(fam, Am_, "float32");
      fwrite(fWo, Wo, "float32");

      if synth_phase

        % synthesis phase spectra from magnitiude spectra using minimum phase techniques

        fft_enc = 512;
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
  
   if exist("voicing_", "var")
     v_out_name = sprintf("%s_v.txt", output_prefix);
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


% Basic unquantised rate K mel-sampling then back to rate L

function [model_ rate_K_surface sd_log] = experiment_mel_freq(model, vq_en=0, plots=1, voicing)
  [frames nc] = size(model);
  K = 20; Fs = 8000;
  AmdB = zeros(frames, 160);
  
  for f=1:frames
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    AmdB(f,1:L) = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*Fs/(2000*pi);
    [rate_K_vec rate_K_sample_freqs_kHz] = resample_const_rate_f_mel(model(f,:), K, Fs, 'lanc');
    rate_K_surface(f,:) = rate_K_vec;
  end

  if plots
    mesh(rate_K_surface);
  end

  [model_ AmdB_ ] = resample_rate_L(model, rate_K_surface, rate_K_sample_freqs_kHz, Fs, 'lancmel');

  % calculate SD

  sd_log = [];
  for f=1:frames
    L = model(f,2);
    asd = std(AmdB(f,1:L) - AmdB_(f,1:L));

    % this code useful to explore outliers, adjust threshold based on plo(sd_log) below
    plot_outliers = 0;
    if plot_outliers && (asd > 7)
      figure; plot(AmdB(f,1:L), 'b+-'); hold on; plot(AmdB_(f,1:L), 'r+-'); hold off;      
    end
    sd_log = [sd_log asd];
  end

  if plots
    figure; plot(sd_log); title('SD againstframe');
    figure; plot(hist(sd_log));
  end
  printf("mean SD %4.2f\n", mean(sd_log));
endfunction


