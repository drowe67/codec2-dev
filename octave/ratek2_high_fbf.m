% ratek2_high_fbf.m
%
% David Rowe 2022
%
% Rate K Experiment 2 - Filtering Am, resampling rate L<->K
%                     - interactive Octave script to explore frame by frame
%                       operation of rate K resampling
%                     - this version upsamples to a rate Lhigh before filtering
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README.md for
%    instructions.
%   ~/codec2-dev/build_linux$ ./c2sim ../raw/big_dog.raw --hpf --dump big_dog
%   $ cd ~/codec2-dev/octave
%   octave:14> ratek2_high_fbf("../build_linux/big_dog",50)


function ratek2_high_fbf(samname, f, vq_stage1_f32="", vq_stage2_f32="")
  more off;

  newamp_700c; melvq;
  Fs = 8000; Nb = 20; K = 30; resampler = 'spline'; Lhigh = 80; vq_en = 0; all_en = 0;

  % load up text files dumped from c2sim ---------------------------------------

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);
  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K);

  % optionally load up VQ

  if length(vq_stage1_f32)
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

  % precompute filters at rate Lhigh. Note range of harmonics is 1:Lhigh-1, as
  % we don't use Lhigh-th harmonic as it's on Fs/2

  h = zeros(Lhigh, Lhigh);
  F0high = (Fs/2)/Lhigh;
  figure(2); clf; hold on;
  for m=1:Lhigh-1
    h(m,:) = generate_filter(m,F0high,Lhigh,Nb);
    plot((1:Lhigh-1)*F0high,h(m,1:Lhigh-1))
  end
  hold off;
  
  rate_Lhigh_sample_freqs_kHz = (F0high:F0high:(Lhigh-1)*F0high)/1000;

  % Keyboard loop --------------------------------------------------------------

  k = ' '; 
  do
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    figure(1); clf; plot(s); axis([1 length(s) -20000 20000]);

    Wo = model(f,1); F0 = Fs*Wo/(2*pi); L = model(f,2);
    Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    % resample from rate L to rate Lhigh (both linearly spaced)
    
    AmdB_rate_Lhigh = interp1([0 Am_freqs_kHz 4], [0 AmdB 0], rate_Lhigh_sample_freqs_kHz, "spline", "extrap");
    
    % Filter at rate Lhigh, y = F(R(a)). Note we filter in linear energy domain, and Lhigh are linearly spaced

    YdB = zeros(1,Lhigh-1);
    for m=1:Lhigh-1
      Am_rate_Lhigh = 10.^(AmdB_rate_Lhigh/20);
      Y = sum(Am_rate_Lhigh.^2 .* h(m,1:Lhigh-1));
      YdB(m) = 10*log10(Y);
    end
    
    % Resample to rate K, optionally VQ, then back to rate Lhigh to check error

    B = interp1(rate_Lhigh_sample_freqs_kHz, YdB, rate_K_sample_freqs_kHz, "spline", "extrap");

    Eq = 0;
    if vq_en
      amean = mean(B);
      %[mse_list index_list] = search_vq(vq, B-amean, 1);
      [res B_hat ind] = mbest(vq, B-amean, mbest_depth);
      B_hat = B_hat + amean;
      Eq = sum((B-B_hat).^2)/K;
      B = B_hat;
    end
    YdB_ = interp1([0 rate_K_sample_freqs_kHz 4], [0 B 0], rate_Lhigh_sample_freqs_kHz, "spline", 0);
    
    figure(3); clf;
    hold on;
    l = sprintf(";rate %d AmdB;g+-", L);
    if all_en
      plot((1:L)*Wo*4000/pi, AmdB, l);
      plot(rate_Lhigh_sample_freqs_kHz*1000, AmdB_rate_Lhigh, ';rate Lhigh AdB;k+-');    
    end
    plot(rate_Lhigh_sample_freqs_kHz*1000, YdB, ';rate Lhigh YdB;b+-');    
    stem(rate_K_sample_freqs_kHz*1000, B, ";rate K;c+-");

    Lmin = round(200/F0high); Lmax = floor(3700/F0high);
    E = sum((YdB(Lmin:Lmax) - YdB_(Lmin:Lmax)).^2)/(Lmax-Lmin+1);
    plot((1:Lhigh-1)*F0high, YdB_,";rate Lhigh YdB hat;r+-");
    l = sprintf(";Eq %3.2f E %3.2f dB;bk+-", Eq, E);
    plot((Lmin:Lmax)*F0high, (YdB(Lmin:Lmax) - YdB_(Lmin:Lmax)), l);
    axis([0 Fs/2 -10 80]);
    hold off;

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit p-png v-vq[%d] a-all plots", f, vq_en);
    fflush(stdout);
    k = kbhit();

    if k == 'n'
      f = f + 1;
    endif
    if k == 'b'
      f = f - 1;
    endif
    if k == 'e'
      if energy == 1
        energy = 0;
      else
        energy = 1;
      end
    end
    if k == 'v'
      if vq_en, vq_en = 0; else vq_en = 1; end
    end
    if k == 'a'
      if all_en, all_en = 0; else all_en = 1; end
    end
    if (k == 'p')
      [dir name ext]=fileparts("../build_linux/big_dog");
      set(gca, 'FontSize', 16);
      h = legend({"Rate L Am","Rate K Bm", "Rate L Am hat"}, "location", "north");
      legend("boxoff")
      set (h, "fontsize", 16);
      xlabel('Freq (Hz)'); ylabel('Amplitude (dB)');
      print(sprintf("ratek2_%s_%d",name,f),"-dpng","-S500,500");
    endif

  until (k == 'q')
  printf("\n");

endfunction
