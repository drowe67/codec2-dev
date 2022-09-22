% ratek2_batch.m
%
% David Rowe 2022
%
% Rate K Experiment 2 - Filtering Am, resampling rate L<->K
%                     - batch mode version
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README.md for
%    instructions.
%   ~/codec2-dev/build_linux$ ./c2sim ../raw/big_dog.raw --dump big_dog
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

  % final curve with para - models current rateK
  K_vec = []; E_vec = [];
  Nb = 100;
  for K=[10 15 20 25 30 40 50]
    E = aratek2_batch(samname, Nb, K, 'para');
    K_vec = [K_vec K];
    E_vec = [E_vec mean(E)];
  end
  plot(K_vec, E_vec,'-+');
  leg{i} = sprintf('Nb=%d para',Nb); i++;
  
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

