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
  ./src/c2sim --hpf ~/raw/big_dog.raw --modelout big_dog
  
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
  
  Nb=20; K=30; rateK_en = 0; verbose = 1; eq =0;
  A_out_fn = ""; B_out_fn = ""; vq_stage1_f32=""; vq_stage2_f32=""; vq_stage3_f32="";
  H_out_fn = ""; amp_pf_en = 0;  phase_pf_en=0; i = 1;
  Kst=0; Ken=K-1; dec = 1; scatter_en = 0; noise_var = 0;
  w = ones(1,K); w1 = ones(1,K); dec_lin = 1; pre_en = 0; logfn="";

  lower = 10;             % only consider vectors above this mean
  dynamic_range = 100;     % restrict dynamic range of vectors
  
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
    elseif strcmp(varargin{i},"vq1") 
      vq_stage1_f32 = varargin{i+1}; i++;
      rateK_en = 1;
    elseif strcmp(varargin{i},"vq2") 
      vq_stage2_f32 = varargin{i+1}; i++;
    elseif strcmp(varargin{i},"vq3") 
      vq_stage3_f32 = varargin{i+1}; i++;
    elseif strcmp(varargin{i},"amp_pf") 
      amp_pf_en = 1;
    elseif strcmp(varargin{i},"phase_pf") 
      phase_pf_en = 1;
    elseif strcmp(varargin{i},"eq") 
      eq = varargin{i+1}; i++;
    elseif strcmp(varargin{i},"verbose") 
      verbose = 2;
    elseif strcmp(varargin{i},"K") 
      rateK_en = 1;
      K = varargin{i+1}; i++;
      Kst=0; Ken=K-1; w = ones(1,K);
      w1 = ones(1,K);
   elseif strcmp(varargin{i},"subset") 
      % restrict range of VQ match to a subset of rate K vector B
      % these are in C index format for compatability with C
      % so default is Kst=0 Ken=K-1;
      Kst = 2; Ken = 24;
      w(1:Kst) = 0; w(Ken+2:K) = 0;
    elseif strcmp(varargin{i},"dec") 
      dec = varargin{i+1}; i++;
    elseif strcmp(varargin{i},"DR") 
      dynamic_range = varargin{i+1}; i++;
    elseif strcmp(varargin{i},"scatter") 
      % energy scatter plots
      scatter_en = 1;
    elseif strcmp(varargin{i},"noise") 
      % add noise to B vector
      noise_var = varargin{i+1}; i++;
    elseif strcmp(varargin{i},"weights") 
      % user defined weights
      w1 = varargin{i+1}; i++;
    elseif strcmp(varargin{i},"nearest") 
      % choose nearest when decimating
      dec_lin = 0;
    elseif strcmp(varargin{i},"pre") 
      pre_en = 1;    
    elseif strcmp(varargin{i},"logfn") 
      logfn = varargin{i+1}; i++;
      printf("logfn: %s\n", logfn);
  else
      printf("\nERROR unknown argument: %s\n", varargin{i});
      return;
    end
    i++;      
  end  

  model_name = strcat(samname,"_model.bin");
  model = load_codec2_model(model_name);
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
    if length(vq_stage3_f32)
      vq_stage3 = load_f32(vq_stage3_f32,K);
      vq(:,:,3)= vq_stage3;
      [M tmp] = size(vq_stage3); printf("stage 3 vq size: %d\n", M);
    end
  end

  sum_Eq = 0; nEq = 0;
  
  for f=1:frames
    Wo = model(f,1); F0 = Fs*Wo/(2*pi); L = model(f,2);
    Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    rate_L_sample_freqs_kHz = ((1:L)*F0)/1000;
    
    % optionally apply pre-emphasis
    if pre_en
      p = 1 - cos(Wo*(1:L)) + j*sin(Wo*(1:L));
      PdB = 20*log10(abs(p));
      AmdB += PdB;
    end
    
    % resample from rate L to rate Lhigh (both linearly spaced)

    AmdB_rate_Lhigh = interp1([0 rate_L_sample_freqs_kHz 4], [0 AmdB 0], rate_Lhigh_sample_freqs_kHz, "spline", "extrap");
    
    % Filter at rate Lhigh, y = F(R(a)). Note we filter in linear energy domain, and Lhigh are linearly spaced

    YdB = zeros(1,Lhigh-1);
    for m=1:Lhigh-1
      Am_rate_Lhigh = 10.^(AmdB_rate_Lhigh/20);
      Y = sum(Am_rate_Lhigh.^2 .* h(m,1:Lhigh-1));
      YdB(m) = 10*log10(Y);
    end
    
    if rateK_en
      % Resample from rate Lhigh to rate K b=R(Y), note K are non-linearly spaced (warped freq axis)
      B(f,:) = interp1(rate_Lhigh_sample_freqs_kHz, YdB, rate_K_sample_freqs_kHz, "spline", "extrap");
      
      % dynamic range limiting
      lower=-100;
      B(f,:) .*= w;
      target = zeros(1,K);
      target(Kst+1:Ken+1) = B(f,Kst+1:Ken+1);
      mx = max(target(Kst+1:Ken+1)); target(Kst+1:Ken+1) = max(target(Kst+1:Ken+1), mx - dynamic_range);
      amean = sum(target)/(Ken-Kst+1);
      target -= amean;
      if vq_en
        [res target_ ind] = mbest(vq, target, mbest_depth, w1);
        Eq(f) = sum((target-target_).^2)/(Ken-Kst+1);
        if amean > lower, sum_Eq += Eq(f); nEq++; end
        if verbose >= 2
          printf("f: %3d amean: %5.1f Eq: %4.2f dB^2", f, amean, Eq(f));
          for i=1:length(ind)
            printf(" %4d",ind(i));
          end
          printf("\n");
        end  
      else
        target_ = target;
      end
      B_hat(f,:) = target_ + amean;
      
      % optional noise injection to simulate quantisation - low freq samples 
      % appear more sensitive to quantisation noise
      B_hat(f,1:K) += sqrt(noise_var)*randn(1,K);

      % ensure frame energy is unchanged after quantisation
      Blin = 10 .^ (B(f,:)/20); E1 = sum(Blin .^2);
      Blin_hat = 10 .^ (B_hat(f,:)/20); E2 = sum(Blin_hat .^2);      
      B_hat(f,:) += 10*log10(E1/E2);
       
      AmdB_ = interp1([0 rate_K_sample_freqs_kHz 4], [0 B_hat(f,:) 0], rate_L_sample_freqs_kHz, "spline", 0);
      nzero = floor(rate_K_sample_freqs_kHz(Kst+1)*1000/F0);
      AmdB_(1:nzero) = 0;
      if Ken+1 < K
        nzero = floor(rate_K_sample_freqs_kHz(Ken+2)*1000/F0);
        AmdB_(nzero:L) = 0;
      end
      
      % Optional amplitude post filtering
      if amp_pf_en
        AmdB_ = amplitude_postfilter(rate_L_sample_freqs_kHz, AmdB_, Fs, F0, eq);
      end

      if pre_en
        AmdB_ -= PdB;
      end
      Am_(f,1:L) = 10.^(AmdB_/20);      
      
      % Synthesised phase0 model using Hilbert Transform
      if length(H_out_fn)
        H(f,1:L) = synth_phase_from_mag(rate_K_sample_freqs_kHz, B_hat(f,:), Fs, Wo, L, phase_pf_en);
      end
    else
      % Optional amplitude post filtering
      if amp_pf_en
        YdB = amplitude_postfilter(rate_Lhigh_sample_freqs_kHz, YdB, Fs, F0high, eq);
      end
      
      % back to rate L
      AmdB_ = interp1([0 rate_Lhigh_sample_freqs_kHz 4], [0 YdB 0], rate_L_sample_freqs_kHz, "spline", "extrap");
      if pre_en
        AmdB_ += PdB;
      end
      Am_(f,1:L) = 10.^(AmdB_/20);
      
      if length(H_out_fn)
        H(f,1:L) = synth_phase_from_mag(rate_Lhigh_sample_freqs_kHz, YdB, Fs, Wo, L, phase_pf_en);
      end
    end
    printf("%d/%d %3.0f%%\r", f,frames, (f/frames)*100);
  end
  printf("\n");
  
  if dec > 1
    for f=1:dec:frames-dec
      fa = f; fb = f+dec;
      x = 1/dec; x_inc = 1/dec;
      for sf=fa+1:fb-1 
        Wo = model(sf,1); F0 = Fs*Wo/(2*pi); L = model(sf,2);
        rate_L_sample_freqs_kHz = ((1:L)*F0)/1000;
        if dec_lin
          % linear interpolation, tends to muffle speech
          B_hat(sf,:) = (1-x)*B_hat(fa,:) + x*B_hat(fb,:);
          printf("f: %d sf: %d fa: %d fb: %d x: %3.2f\n", f,sf, fa,fb,x);
        else
          % choose closest, requires an extra bit, doesn't sound very good (could be a bug)
          dista = sum((B_hat(sf,:) - B_hat(fa,:)) .^ 2);
          distb = sum((B_hat(sf,:) - B_hat(fb,:)) .^ 2);
          if dista < distb
            B_hat(sf,:) = B_hat(fa,:);
            choice = 'a';
          else
            B_hat(sf,:) = B_hat(fb,:);
            choice = 'b';
          end
          printf("f: %d sf: %d fa: %d fb: %d choice: %c\n", f, sf, fa, fb, choice);
        end
        
        x += x_inc;
        AmdB_ = interp1([0 rate_K_sample_freqs_kHz 4], [0 B_hat(sf,:) 0], rate_L_sample_freqs_kHz, "spline", 0);
        
        % Optional amplitude post filtering
        if amp_pf_en
          AmdB_ = amplitude_postfilter(rate_L_sample_freqs_kHz, AmdB_, Fs, F0, eq);
        end
        nzero = floor(rate_K_sample_freqs_kHz(Kst+1)*1000/F0);
        AmdB_(1:nzero) = 0;
        if Ken+1 < K
          nzero = floor(rate_K_sample_freqs_kHz(Ken+2)*1000/F0);
          AmdB_(nzero:L) = 0;
        end
        Am_(sf,1:L) = 10.^(AmdB_/20);
 
        if length(H_out_fn)
          H(sf,1:L) = synth_phase_from_mag(rate_K_sample_freqs_kHz, B_hat(sf,:), Fs, Wo, L, phase_pf_en);
        end
      end
    end
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
    sd = sum_Eq/nEq;
    printf("Nb: %d K: %d mean Eq: %4.2f dB^2 %4.2f dB\n", Nb, K, sd, sqrt(sd));
    if length(logfn)
      fl = fopen(logfn,"wt"); fprintf(fl,"%4.2f\n",sd); fclose(fl);
    end
  end
  
  % optional energy scatter plots
  if scatter_en
    E_A = zeros(1,frames); E_B = zeros(1,frames); Wo = zeros(1,frames);
    for f=1:frames
      Wo(f) = model(f,1); L = model(f,2);
      Am = model(f,3:(L+2));
      E_A(f) = 10*log10(sum(Am.^2));
      Blin = 10 .^ (B(f,:)/20);
      % try compensating for higher energy high F0 Bm
      %Blin *= L/K;
      E_B(f) = 10*log10(sum(Blin .^2));
    end
    figure(1); clf; 
    semilogx(Fs*Wo/(2*pi),E_A,'b.;Energy {A};');
    hold on;
    semilogx(Fs*Wo/(2*pi),E_B,'g.;Energy {B};');
    hold off; axis([50 400 20 80])
  end
endfunction

