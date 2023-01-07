% ratek2_high_fbf.m
%
% David Rowe 2022
%
% Rate K Experiment 2 - Resample rate L to rate Lhigh, filter, and optionally VQ
%                     - interactive Octave script to explore frame by frame
%                       operation of rate K resampling
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
  amp_pf_en = 1; eq = 0; Kst = 2; Ken = 24;
  
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
    vq_en = 1;
    vq_stage1 = load_f32(vq_stage1_f32,K);
    vq(:,:,1)= vq_stage1; 
    [M tmp] = size(vq_stage1); printf("stage 1 vq size: %d\n", M);
    nvq = 1; mbest_depth = 1;
    if length(vq_stage2_f32)
      vq_stage2 = load_f32(vq_stage2_f32,K);
      vq(:,:,2)= vq_stage2; 
      [M tmp] = size(vq_stage2); printf("stage 2 vq size: %d\n", M);
      nvq++; mbest_depth = 5;
    end
    w = zeros(1,K); w(Kst+1:Ken+1) = 1
  end

  % precompute filters at rate Lhigh. Note range of harmonics is 1:Lhigh-1, as
  % we don't use Lhigh-th harmonic as it's on Fs/2

  h = zeros(Lhigh, Lhigh);
  F0high = (Fs/2)/Lhigh;
  for m=1:Lhigh-1
    h(m,:) = generate_filter(m,F0high,Lhigh,Nb);
    plot((1:Lhigh-1)*F0high,h(m,1:Lhigh-1))
  end
  
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
    
    % Optional amplitude post filtering
    if amp_pf_en
      [YdB SdB] = amplitude_postfilter(rate_Lhigh_sample_freqs_kHz, YdB, Fs, F0high, eq);
      if vq_en == 0
        figure(2); clf;
        plot(rate_Lhigh_sample_freqs_kHz*1000, YdB-SdB);
        axis([0 4000 -20 20]);
      end        
    end
   
    % Resample to rate K, optionally VQ, then back to rate Lhigh to check error
    B = interp1(rate_Lhigh_sample_freqs_kHz, YdB, rate_K_sample_freqs_kHz, "spline", "extrap");
    
    Eq = 0;
    if vq_en
      #m = max(B); B = max(B, m-40);
      B .*= w;
      amean = sum(B)/(Ken-Kst+1);
      target = zeros(1,K);
      target(Kst+1:Ken+1) = B(Kst+1:Ken+1)-amean;
      [res target_ ind] = mbest(vq, target, mbest_depth);
      B_hat = target_; B_hat(Kst+1:Ken+1) += amean;
      Eq = sum((target-target_).^2)/(Ken-Kst+1);
      figure(2); clf; hold on;
      for i=1:nvq
        plot(rate_K_sample_freqs_kHz*1000, vq(ind(i),:,i),'b;vq;');
      end
      plot(rate_K_sample_freqs_kHz*1000, target,'g;target;');
      B = B_hat
      plot([0 4000], [amean amean]);
      hold off; axis([0 4000 -40 60]);
    end
    YdB_ = interp1([0 rate_K_sample_freqs_kHz 4], [0 B 0], rate_Lhigh_sample_freqs_kHz, "spline", 0);
    nzero = floor(rate_K_sample_freqs_kHz(Kst)*1000/F0high);
    YdB_(1:nzero) = 0;
    
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
    le = sprintf("Eq %3.2f E %3.2f dB", Eq, E);
    plot((Lmin:Lmax)*F0high, (YdB(Lmin:Lmax) - YdB_(Lmin:Lmax)), sprintf(";%s;bk+-",le));
    axis([0 Fs/2 -10 80]);
    hold off;

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit p-png f-postfilter e-eq[%d] v-vq[%d] a-all plots", f, eq, vq_en);
    if vq_en 
      printf(" Eq: %3.2f dB^2", Eq);
      for i=1:length(ind)
        printf(" %4d",ind(i));
      end
    end  
    fflush(stdout);
    k = kbhit();

    if k == 'n', f = f + 1; end
    if k == 'b', f = f - 1; end
    if k == 'v', vq_en = mod(vq_en+1,2); end
    if k == 'a', all_en = mod(all_en+1,2); end
    if k == 'e', eq = mod(eq+1,3); end
    if k == 'f', amp_pf_en = mod(amp_pf_en+1,2); end
    if (k == 'p')
      [dir name ext]=fileparts(samname);
      set(gca, 'FontSize', 16);
      hl = legend({"Rate Lhigh YdB","Rate K" "Rate Lhigh YdB hat",le}, "location", "northeast");
      legend("boxoff")
      set (hl, "fontsize", 16);
      xlabel('Freq (Hz)'); ylabel('Amplitude (dB)');
      print(sprintf("ratek2_%s_%d",name,f),"-dpng","-S500,500");
    endif

  until (k == 'q')
  printf("\n");

endfunction
