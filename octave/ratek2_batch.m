% ratek2_batch.m
%
% David Rowe 2022
%
% Rate K Experiment 2 - Filtering Am, resampling rate L<->K
%                     - batch mode version, various functions
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README.md for
%    instructions.
%   ~/codec2-dev/build_linux$ ./src/c2sim ../raw/big_dog.raw --dump big_dog
%   $ cd ~/codec2-dev/octave
%   octave:14> ratek2_batch("../build_linux/big_dog")

1;

function ratek2_batch_curves(samname)
  more off;

  h=figure(1); clf; hold on; i=1;

  for Nb=[10 15 20 25 100]
    K_vec = []; E_vec = [];
    for K=[10 15 20 25 30 40 50]
      E = aratek2_batch(samname, Nb, K);
      K_vec = [K_vec K];
      E_vec = [E_vec mean(E)];
    end
    plot(K_vec, E_vec,'-+');
    leg{i} = sprintf('Nb=%d',Nb); i++;
  end

  % curve with para - models current rateK
  K_vec = []; E_vec = [];
  Nb = 100;
  for K=[10 15 20 25 30 40 50]
    E = aratek2_batch(samname, Nb, K, 'para');
    K_vec = [K_vec K];
    E_vec = [E_vec mean(E)];
  end
  plot(K_vec, E_vec,'-+');
  leg{i} = sprintf('Nb=%d para',Nb); i++;

  % curve with efficient ratek2
  K_vec = []; E_vec = [];
  Nb = 20;
  for K=[10 15 20 25 30 40 50]
    E = aratek2_batch_efficient(samname, Nb, K);
    K_vec = [K_vec K];
    E_vec = [E_vec mean(E)];
  end
  plot(K_vec, E_vec,'-o');
  leg{i} = sprintf('Nb=%d Eff',Nb); i++;

  hold off;
  axis([0 50 0 10]);
  xlabel('K'); ylabel('mean E (dB^2)');
  set(gca, 'FontSize', 16);
  h = legend(leg, "location", "east");
  legend("boxoff");
  set (h, "fontsize", 16);
  [dir name ext] = fileparts(samname);
  print("-dpng", sprintf("ratek2_E_K_%s",name), "-S500,500");

  % TODO E against F0 (any F0 depedancy)
  % TODO E against frame (and big outliers)

endfunction


function [E F0] = aratek2_batch(samname, Nb=20, K=30, resampler='spline', Y_out_fn="")
  more off;

  newamp_700c;
  Fs = 8000; max_amp = 160;

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K);
  E = zeros(1,frames);
  F0 = zeros(1,frames);
  Y = zeros(frames, max_amp);

  for f=1:frames
    Wo = model(f,1); F0(f) = Fs*Wo/(2*pi); L = model(f,2);
    Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    % Filter at rate L, Y = F(A)
    Y = zeros(1,L); YdB = zeros(1,L);
    for m=1:L
      h = generate_filter(m,F0(f),L,Nb);
      Y(m) = sum(Am.^2 .* h);
      YdB(m) = 10*log10(Y(m));
    end

    % Resample to rate K, then back to rate L
    amodel = model(f,:); amodel(3:(L+2)) = 10 .^ (YdB/20);
    B = resample_const_rate_f(amodel, rate_K_sample_freqs_kHz, clip_en = 0, resampler);
    model_ = resample_rate_L(amodel, B, rate_K_sample_freqs_kHz, resampler);
    Y_(f,1:L) = model_(1,3:(L+2)); YdB_ = 20*log10(Y_(f,1:L));

    Lmin = round(200/F0(f)); Lmax = floor(3700/F0(f));
    E(f)  = sum((YdB(Lmin:Lmax) - YdB_(Lmin:Lmax)).^2)/(Lmax-Lmin+1);
  end

  % optionally write Y, so we can listen using:
  %   ./src/c2sim ../raw/big_dog.raw --amread ../octave/big_dog_y.f32 -o - | aplay -f S16_LE
  if length(Y_out_fn)
    max_amp = 160;
    fy = fopen(Y_out_fn,"wb");
    for f=1:frames
      Yfloat_ = zeros(1,max_amp);
      L = model(f,2);
      Yfloat_(2:L+1) = Y_(f,1:L);
      fwrite(fy, Yfloat_, "float32");
    end
    fclose(fy);
  end

  printf("Nb: %d K: %d mean SD: %4.2f dB^2\n", Nb, K, mean(E));
endfunction

% More efficient implementation, we resample from rate L to Lhigh then filter at rate Lhigh,
% so that we only need to calculate filters once.

function [E F0] = aratek2_batch_efficient(samname, Nb=20, K=30)
  more off;

  newamp_700c;
  Fs = 8000; max_amp = 160;  Lmax = 80; Lhigh = Lmax;

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames model_cols] = size(model);
  amodel = zeros(frames,model_cols);
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K);
  E = zeros(1,frames);
  F0 = zeros(1,frames);

  % precompute filters at rate Lhigh. Note range of harmonics is 1:Lhigh-1, as
  % we don't use Lhigh-th harmonic as it's on Fs/2

  h = zeros(Lhigh, Lhigh);
  F0high = (Fs/2)/Lhigh;
  for m=1:Lhigh-1
    h(m,:) = generate_filter(m,F0high,Lhigh,Nb);
  end
  rate_Lhigh_sample_freqs_kHz = (F0high:F0high:(Lhigh-1)*F0high)/1000;

  for f=1:frames
    Wo = model(f,1); F0(f) = Fs*Wo/(2*pi); L = model(f,2);
    Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    rate_L_sample_freqs_kHz = ((1:L)*F0(f))/1000;

    % resample from rate L to rate Lhigh (both linearly spaced)

    AmdB_rate_Lhigh = interp1([0 rate_L_sample_freqs_kHz 4], [0 AmdB 0], rate_Lhigh_sample_freqs_kHz, "spline", "extrap");

    % Filter at rate Lhigh, y = F(R(a)). Note we filter in linear energy domain, and Lhigh are linearly spaced

    YdB = zeros(1,Lhigh-1);
    for m=1:Lhigh-1
      Am_rate_Lhigh = 10.^(AmdB_rate_Lhigh/20);
      Y = sum(Am_rate_Lhigh.^2 .* h(m,1:Lhigh-1));
      YdB(m) = 10*log10(Y);
    end

    % Resample from rate Lhigh to rate K b=R(Y), note K are non-linearly spaced (warped freq axis)

    B = interp1(rate_Lhigh_sample_freqs_kHz, YdB, rate_K_sample_freqs_kHz, "spline", "extrap");

    % now back to rate Lhigh

    YdB_ = interp1([0 rate_K_sample_freqs_kHz 4], [0 B 0], rate_Lhigh_sample_freqs_kHz, "spline", 0);

    Lmin = round(200/F0high); Lmax = floor(3700/F0high);
    E(f)  = sum((YdB(Lmin:Lmax) - YdB_(Lmin:Lmax)).^2)/(Lmax-Lmin+1);
  end

  printf("Nb: %d K: %d mean SD: %4.2f dB^2\n", Nb, K, mean(E));
endfunction


% used to generate rate K VQ training data
#{
  To compare VQ perf with Nb=20 and Nb=100:

  cd codec2/build_linux
  ./src/c2sim ~/Downloads/train_120.spc --dump train_120
  octave:49> ratek2_batch;  B = ratek2_model_to_ratek("../build_linux/train_120",20,30,"train_120_Nb20_K30.f32");
  octave:50> ratek2_batch;  B = ratek2_model_to_ratek("../build_linux/train_120",100,30,"train_120_Nb100_K30.f32");
  $ ../script/ratek_resampler.sh ../octave/train_120_Nb20_K30.f32
  $ ../script/ratek_resampler.sh ../octave/train_120_Nb100_K30.f32
#}

function B = ratek2_model_to_ratek(samname, Nb=20, K=30, B_out_fn="", vq_stage1_f32="", vq_stage2_f32="", Y_out_fn="")
  more off;

  newamp_700c;
  Fs = 8000; max_amp = 160; resampler='spline'; Lhigh=80; max_amp = 160;

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

    % Resample from rate Lhigh to rate K b=R(Y), note K are non-linearly spaced (warped freq axis)
    B(f,:) = interp1(rate_Lhigh_sample_freqs_kHz, YdB, rate_K_sample_freqs_kHz, "spline", "extrap");

    if vq_en
      amean = mean(B(f,:));
      [res B_hat ind] = mbest(vq, B(f,:)-amean, mbest_depth);
      B_hat = B_hat + amean;
      Eq(f) = sum((B(f,:)-B_hat).^2)/K;
      YdB_ = interp1([0 rate_K_sample_freqs_kHz 4], [0 B_hat 0], rate_L_sample_freqs_kHz, "spline", 0);
      Y_(f,1:L) = 10.^(YdB_/20);
      printf("f: %d Eq: %3.2f dB\n", f, Eq(f));
    else
      YdB_ = interp1([0 rate_K_sample_freqs_kHz 4], [0 B(f,:) 0], rate_L_sample_freqs_kHz, "spline", 0);
      Y_(f,1:L) = 10.^(YdB_/20);
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

  % optionally write Y, so we can listen using:
  %   ./src/c2sim ../raw/big_dog.raw --amread ../octave/big_dog_y.f32 -o - | aplay -f S16_LE
  if length(Y_out_fn)
     fy = fopen(Y_out_fn,"wb");
    for f=1:frames
      Yfloat_ = zeros(1,max_amp);
      L = model(f,2);
      Yfloat_(2:L+1) = Y_(f,1:L);
      fwrite(fy, Yfloat_, "float32");
    end
    fclose(fy);
  end

  if vq_en
    printf("Nb: %d K: %d mean SD: %4.2f dB^2\n", Nb, K, mean(Eq));
  end
endfunction

#{
  Experiment to test post filter algorithm developed in plphase2.m (fbf version).
  At rate Lhigh, Nb filter, then apply amplitude and phase post filters.

  ./src/c2sim ../raw/big_dog.raw --dump big_dog
  octave:90> ratek2_batch;  ratek2_model_postfilter("../build_linux/big_dog","big_dog_a.f32","big_dog_h.f32",1);
  ./src/c2sim ../raw/big_dog.raw --phase0 --postfilter --amread ../octave/big_dog_a.f32 -o - | aplay -f S16_LE
  ./src/c2sim ../raw/big_dog.raw --phase0 --postfilter --amread ../octave/big_dog_a.f32 --hmread ../octave/big_dog_h.f32 -o - | aplay -f S16_LE
#}

function ratek2_model_postfilter(samname, A_out_fn="", H_out_fn="", amp_pf_en=0, phase_pf_en=0)
  more off;

  newamp_700c;
  Fs = 8000; max_amp = 160; resampler='spline'; Lhigh=80; max_amp = 160;
  Nb=20; K=30;

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K);

  % precompute filters at rate Lhigh. Note range of harmonics is 1:Lhigh-1, as
  % we don't use Lhigh-th harmonic as it's on Fs/2

  h = zeros(Lhigh, Lhigh);
  F0high = (Fs/2)/Lhigh;
  for m=1:Lhigh-1
    h(m,:) = generate_filter(m,F0high,Lhigh,Nb);
  end
  rate_Lhigh_sample_freqs_kHz = (F0high:F0high:(Lhigh-1)*F0high)/1000;

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

    if amp_pf_en
      YdB = amplitude_postfilter(rate_Lhigh_sample_freqs_kHz, YdB, Fs, F0high);
    end
    % Synthesised phase0 model using Hilbert Transform
    H(f,1:L) = synth_phase_from_mag(rate_Lhigh_sample_freqs_kHz, YdB, Fs, Wo, L, phase_pf_en);

    AmdB_ = interp1([0 rate_Lhigh_sample_freqs_kHz 4], [0 YdB 0], rate_L_sample_freqs_kHz, "spline", "extrap");
    Am_(f,1:L) = 10.^(AmdB_/20);
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
endfunction
