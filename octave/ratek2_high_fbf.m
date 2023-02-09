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
  Fs = 8000; Nb = 20; K = 20; resampler = 'spline'; Lhigh = 80; vq_en = 0; all_en = 0;
  amp_pf_en = 0; eq = 0; Kst = 0; Ken = K-1; pre_en = 0; w_en = 0;
  
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
    w = ones(1,K);
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
    
    % optionally apply pre-emphasis
    if pre_en
      p = 1 - cos(Wo*(1:L)) + j*sin(Wo*(1:L));
      PdB = 20*log10(abs(p));
      AmdB += PdB;
    end
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
      if w_en
        % weighted search, requires gain calculation for each 
        % codebook entry. we only support single stage
        mx = max(target);
        w = (0.75/30)*(target-mx) + 1.0;
        w2 = w .^ 2;
        I = length(vq_stage1);
        Emin = 1E32;
        for i=1:I
          g = sum((B-vq_stage1(i,:)).*w2)/sum(w2);
          E = sum( ((B - vq_stage1(i,:) - g) .* w) .^ 2 );
          if E < Emin
            best_i = i;
            Emin = E;
            gmin = g;
          end
        end
        Eq = Emin/K; B_hat = vq_stage1(best_i,:) + gmin;

        figure(2); clf; hold on;
        plot(rate_K_sample_freqs_kHz*1000, vq_stage1(best_i,:),'b;vq;');
        plot(rate_K_sample_freqs_kHz*1000, B-gmin,'g;B-gmin;');
        plot([0 4000], [gmin gmin],'r;gmin;');
        hold off; axis([0 4000 -40 60]);
      else
        % regular unweighted search, we can remove gain/mean outside of loop
        amean = sum(B(Kst+1:Ken+1))/(Ken-Kst+1);
        target = zeros(1,K);
        target(Kst+1:Ken+1) = B(Kst+1:Ken+1)-amean;
        [res target_ ind] = mbest(vq, target, mbest_depth);
        B_hat = target_; B_hat(Kst+1:Ken+1) += amean;
        Eq = sum((target-target_).^2)/(Ken-Kst+1);
        gmin = amean; best_i = ind(1);

        figure(2); clf; hold on;
        for i=1:nvq
          plot(rate_K_sample_freqs_kHz*1000, vq(ind(i),:,i),'b;vq;');
        end
        plot(rate_K_sample_freqs_kHz*1000, target,'g;target;');
        plot([0 4000], [amean amean],'r;mean;');
        hold off; axis([0 4000 -40 60]);
      end
      
      B = B_hat;
    end
    
    YdB_ = interp1([0 rate_K_sample_freqs_kHz 4], [0 B 0], rate_Lhigh_sample_freqs_kHz, "spline", 0);
    nzero = floor(rate_K_sample_freqs_kHz(Kst+1)*1000/F0high);
    YdB_(1:nzero) = 0;
    
    % Optional amplitude post filtering
    if amp_pf_en
      [YdB_ SdB] = amplitude_postfilter(rate_Lhigh_sample_freqs_kHz, YdB_, Fs, F0high, eq);
    end  
   
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

    printf("\rframe: %d  menu: n-nxt b-bck q-qt p-pre f-pf w[%d]-wght v-vq[%d] a-all", f, w_en, vq_en);
    if vq_en
      printf(" Eq: %5.2f gmin: %3.1f best_i: %d\n", Eq, gmin, best_i);
    end
    fflush(stdout);
    k = kbhit();

    if k == 'n', f = f + 1; end
    if k == 'b', f = f - 1; end
    if k == 'v', vq_en = mod(vq_en+1,2); end
    if k == 'a', all_en = mod(all_en+1,2); end
    if k == 'e', eq = mod(eq+1,3); end
    if k == 'f', amp_pf_en = mod(amp_pf_en+1,2); end
    if k == 'p', pre_en = mod(pre_en+1,2); end
    if k == 'w', w_en = mod(w_en+1,2); end

  until (k == 'q')
  printf("\n");

endfunction
