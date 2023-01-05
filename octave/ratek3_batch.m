% ratek3_batch.m
%
% David Rowe 2022
%
% Rate K Experiment 3 - batch processing tool
%                     - similar signal processing to pl2phase.m

#{
  Ratek tool with various options to support
  + filtering 
  + resampling
  + generation of VQ training data
  + amplitude and phase post filtering
  + files of vectors that can be input to c2sim for synthesis

  cd ~/codec2/build_linux
  ./src/c2sim --hpf ~/raw/big_dog.raw --dump big_dog
  
  cd ~/codec2/octave
  octave:49> ratek3_batch; 
  octave:50> B=ratek3_batch_tool("../build_linux/big_dog","B_out","big_dog_B.f32");
  octave:51> mesh(B);
  
  octave:8> ratek3_batch; ratek3_batch_tool("../build_linux/big_dog","A_out","big_dog_a.f32","H_out","big_dog_h.f32");
  ./src/c2sim --hpf ../raw/big_dog.raw --phase0 --postfilter --amread ../octave/big_dog_a.f32 --hmread ../octave/big_dog_h.f32 -o - | aplay -f S16_LE
  
#}

1;

function B = ratek3_batch_tool(samname, varargin)
  more off;

  newamp_700c;
  Fs = 8000; max_amp = 160; resampler='spline'; Lhigh=80; max_amp = 160;
  
  Nb=20; K=30; rateK_en = 0; verbose = 1; restore_slope = 1;
  A_out_fn = ""; B_out_fn = ""; vq_stage1_f32=""; vq_stage2_f32="";
  H_out_fn = ""; amp_pf_en = 0;  phase_pf_en=0; i = 1;
  while i<=length(varargin)
    if strcmp(varargin{i},"rateK") 
      rateK_en = 1;
    elseif strcmp(varargin{i},"A_out") 
      A_out_fn = varargin{i+1}; i++;
    elseif strcmp(varargin{i},"B_out")
      rateK_en = 1;
      B_out_fn = varargin{i+1}; i++;
    elseif strcmp(varargin{i},"H_out")
      H_out_fn = varargin{i+1}; i++;
    elseif strcmp(varargin{i},"vq_stage1") 
      vq_stage1_f32 = varargin{i+1}; i++;
    elseif strcmp(varargin{i},"vq_stage2") 
      vq_stage2_f32 = varargin{i+1}; i++;
    elseif strcmp(varargin{i},"amp_pf") 
      amp_pf_en = 1;
    elseif strcmp(varargin{i},"phase_pf") 
      phase_pf_en = 1;
    elseif strcmp(varargin{i},"eq") 
      restore_slope = 0;
    else
      printf("\nERROR unknown argument: %s\n", varargin{i});
      return;
    end
    i++;      
  end  

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K);
  B = zeros(frames,K);

  % precompute filters at rate Lhigh. Note range of harmonics is 1:Lhigh-1, as
  % we don't use Lhigh-th harmonic as it's on Fs/2

  h = zeros(Lhigh, Lhigh);
  F0high = (Fs/2)/Lhigh;
  for m=1:Lhigh-1
    h(m,:) = generate_filter(m,F0high,Lhigh,Nb);
  end
  rate_Lhigh_sample_freqs_kHz = (F0high:F0high:(Lhigh-1)*F0high)/1000;

  % optionally load up VQ

  vq_en = 0;
  if length(vq_stage1_f32)
    vq_en = 1;
    vq_stage1 = load_f32(vq_stage1_f32,K);
    vq(:,:,1)= vq_stage1;
    [M tmp] = size(vq_stage1); printf("stage 1 vq size: %d\n", M);
    mbest_depth = 1;
    if length(vq_stage2_f32)
      vq_stage2 = load_f32(vq_stage2_f32,K);
      vq(:,:,2)= vq_stage2;
      [M tmp] = size(vq_stage2); printf("stage 2 vq size: %d\n", M);
      mbest_depth = 5;
    end
  end

  for f=1:frames
    Wo = model(f,1); F0 = Fs*Wo/(2*pi); L = model(f,2);
    Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    rate_L_sample_freqs_kHz = ((1:L)*F0)/1000;

    % resample from rate L to rate Lhigh (both linearly spaced)

    AmdB_rate_Lhigh = interp1([0 rate_L_sample_freqs_kHz 4], [0 AmdB 0], rate_Lhigh_sample_freqs_kHz, "spline", "extrap");

    % Filter at rate Lhigh, y = F(R(a)). Note we filter in linear energy domain, and Lhigh are linearly spaced

    YdB = zeros(1,Lhigh-1);
    for m=1:Lhigh-1
      Am_rate_Lhigh = 10.^(AmdB_rate_Lhigh/20);
      Y = sum(Am_rate_Lhigh.^2 .* h(m,1:Lhigh-1));
      YdB(m) = 10*log10(Y);
    end

    % Optional amplitude post filtering
    if amp_pf_en
      YdB = amplitude_postfilter(rate_Lhigh_sample_freqs_kHz, YdB, Fs, F0high, restore_slope);
    end
    
    if rateK_en
      % Resample from rate Lhigh to rate K b=R(Y), note K are non-linearly spaced (warped freq axis)
      B(f,:) = interp1(rate_Lhigh_sample_freqs_kHz, YdB, rate_K_sample_freqs_kHz, "spline", "extrap");

      if vq_en
        amean = mean(B(f,:));
        [res aB_hat ind] = mbest(vq, B(f,:)-amean, mbest_depth);
        B_hat(f,:) = aB_hat + amean;
        Eq(f) = sum((B(f,:)-B_hat(f,:)).^2)/K;
        YdB_ = interp1([0 rate_K_sample_freqs_kHz 4], [0 B_hat(f,:) 0], rate_L_sample_freqs_kHz, "spline", 0);
        Y_(f,1:L) = 10.^(YdB_/20);
        if verbose >= 2
          printf("f: %d Eq: %3.2f dB\n", f, Eq(f));
        end  
      else
        B_hat(f,:) = B(f,:);
      end
    
      AmdB_ = interp1([0 rate_K_sample_freqs_kHz 4], [0 B_hat(f,:) 0], rate_L_sample_freqs_kHz, "spline", 0);
      Am_(f,1:L) = 10.^(AmdB_/20);
      
      % Synthesised phase0 model using Hilbert Transform
      if length(H_out_fn)
        H(f,1:L) = synth_phase_from_mag(rate_K_sample_freqs_kHz, B_hat(f,:), Fs, Wo, L, phase_pf_en);
      end
    else
      % rate Lhigh processing
      AmdB_ = interp1([0 rate_Lhigh_sample_freqs_kHz 4], [0 YdB 0], rate_L_sample_freqs_kHz, "spline", "extrap");
      Am_(f,1:L) = 10.^(AmdB_/20);
      if length(H_out_fn)
        H(f,1:L) = synth_phase_from_mag(rate_Lhigh_sample_freqs_kHz, YdB, Fs, Wo, L, phase_pf_en);
      end
    end
    printf("%d/%d %3.0f%%\r", f,frames, (f/frames)*100);
  end

  % optionally write B to a .f32 file for external VQ training
  if length(B_out_fn)
    fb = fopen(B_out_fn,"wb");
    for f=1:frames
      Bfloat = B(f,:);
      fwrite(fb, Bfloat, "float32");
    end
    fclose(fb);
  end
  
  % optionally write Am, so we can listen using:
  %   ./src/c2sim ../raw/big_dog.raw --amread ../octave/big_dog_am.f32 -o - | aplay -f S16_LE
  if length(A_out_fn)
    fam = fopen(A_out_fn,"wb");
    for f=1:frames
      Afloat_ = zeros(1,max_amp);
      L = model(f,2);
      Afloat_(2:L+1) = Am_(f,1:L);
      fwrite(fam, Afloat_, "float32");
    end
    fclose(fam);
  end
  
  if length(H_out_fn)
    max_amp = 160;
    fhm = fopen(H_out_fn,"wb");
    for f=1:frames
      Hfloat = zeros(1,2*max_amp);
      L = model(f,2);
      for m=1:L
        Hfloat(2*m+1) = cos(H(f,m));
        Hfloat(2*m+2) = sin(H(f,m));
      end
      fwrite(fhm, Hfloat, "float32");
    end
    fclose(fhm);
  end
  
  if vq_en && verbose
    sd = mean(Eq);
    printf("Nb: %d K: %d mean SD: %4.2f dB^2 %4.2f dB\n", Nb, K, sd, sqrt(sd));
  end
endfunction
