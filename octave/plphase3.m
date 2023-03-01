% plphase3.m
%
% fbf script to explore issues with clicks in synthesied audio.

#{
  Usage:

    $ cd codec2/build_linux
    $ ./src/c2sim ../raw/big_dog.raw --hpf --phase0 --dump big_dog

    $ cd codec2/build_linux/octave
    $ octave-cli
    octave:> plphase3("../build_linux/two_lines", 245)
#}

function plphase3(samname, f, Nb=20, K=20)
  [dir basename ext] = fileparts(samname);

  newamp_700c; melvq;
  Fs = 8000; Fs2 = Fs/2; resampler = 'spline'; Lhigh = 80;

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);
  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  phase_name = strcat(samname,"_phase.txt");
  if (file_in_path(".",phase_name))
    phase = unwrap(load(phase_name),pi,2);
  endif
  Wo_vec = model(:,1); F0_vec = Fs*Wo_vec/(2*pi);
  snr_name = strcat(samname,"_snr.txt");
  snr = load(snr_name);
  
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

  k = ' '; plot_synth_sn=1; phase0_en=1; postfilter_en = 1; ratek_en = 1;
  rand_50 = 0;
  do
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    figure(1); clf; plot(s); axis([1 length(s) -20000 20000]);
    
    Wo = model(f,1); F0 = Fs*Wo/(2*pi); L = model(f,2);
    Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    % resample from rate L to rate Lhigh (both linearly spaced)
    AmdB_rate_Lhigh = interp1([0 Am_freqs_kHz 4], [0 AmdB 0],
                              rate_Lhigh_sample_freqs_kHz, "spline", "extrap");

    if ratek_en
      [YdB Y] = filter_rate_Lhigh(Lhigh,h,AmdB_rate_Lhigh);
    else
      YdB = AmdB_rate_Lhigh; Y = 10 .^ (YdB/20);
    end

    if postfilter_en
      YdB_orig = YdB;
      [YdB SdB] = amplitude_postfilter(rate_Lhigh_sample_freqs_kHz, YdB, Fs, F0high);
    end

    % Synthesised phase0 model using Hilbert Transform
    phase0(f,1:L) = synth_phase_from_mag(rate_Lhigh_sample_freqs_kHz, YdB, Fs, Wo, L, postfilter_en);

    % resample from rate Lhigh to rate L (both linearly spaced)
    AmdB_ = interp1([0 rate_Lhigh_sample_freqs_kHz 4], [0 YdB 0], Am_freqs_kHz, "spline", "extrap");
    Am_ = 10 .^ (AmdB_/20);

    % plot time domain speech ---------------------------------------------
    figure(3); clf;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    if phase0_en, phase_ratek(f,1:L) = phase0(f,1:L); else phase_ratek(f,1:L) = phase(f,1:L); end
    if plot_synth_sn
       N=320;
       if rand_50 && F0 < 60
         printf("rand phases...\n");
         phase_ratek(f,round(L/4):L) += 2*pi*rand(1,L-round(L/4)+1);
       end
       s = synth_time(Wo, L, Am, phase_ratek(f,:), N);
       plot(real(s),sprintf('g;%s;',papr(s)));
    end

    figure(2); clf;
    if postfilter_en == 0
      plot((1:L)*Wo*4000/pi, 20*log10(Am),"g+-;Am;");
    endif
    hold on;
    if postfilter_en
      plot(rate_Lhigh_sample_freqs_kHz*1000, YdB_orig, ';rate Lhigh YdB;b-');
      plot(rate_Lhigh_sample_freqs_kHz*1000, YdB, ';rate Lhigh YdB postfilter;g-');
      plot(rate_Lhigh_sample_freqs_kHz*1000, SdB, ';rate Lhigh SdB;g--');
      plot(rate_Lhigh_sample_freqs_kHz*1000, YdB_orig-SdB, ';rate Lhigh YdB-SdB;-');
    else
      plot(rate_Lhigh_sample_freqs_kHz*1000, YdB, ';rate Lhigh YdB;b-');
    end

    figure(4); clf;
    subplot(211); plot(-2:2,F0_vec(f-2:f+2),sprintf('b+-;F0: %2.0f +/- 2 frames;',F0));
    axis([-2 2 50 400]);
    subplot(212); plot(-2:2,snr(f-2:f+2),'b+-;snr +/- 2 frames;');
    axis([-2 2 0 20]);
    
    % interactive menu

    if phase0_en; s2="orig/[phase0]"; else s2="[orig]/phase0"; end
    printf("\rframe: %d  n-nxt b-bk ratek-r 0-%s f-pf[%d] 5-rand50[%d] q-quit",
           f,s2,postfilter_en, rand_50);
    fflush(stdout);
    k = kbhit();
    if (k == 'n')
      f = f + 1;
    endif
    if (k == 'b')
      f = f - 1;
    endif
    if k == '0',
      if phase0_en, phase0_en = 0; else phase0_en = 1; end
    end
    if k == 'f',
      if postfilter_en, postfilter_en = 0; else postfilter_en = 1; end
    end
    if k == 'r',
      if ratek_en, ratek_en = 0; else ratek_en = 1; end
    end
    if k == '5',
      if rand_50, rand_50 = 0; else rand_50 = 1; end
    end
    until (k == 'q')
  printf("\n");

endfunction

% Filter at rate Lhigh, y = F(R(a)). Note we filter in linear energy domain,
% and Lhigh are linearly spaced
function [YdB Y] = filter_rate_Lhigh(Lhigh, h, AmdB_rate_Lhigh)
  Y = zeros(1,Lhigh-1); YdB = zeros(1,Lhigh-1);
  for m=1:Lhigh-1
    Am_rate_Lhigh = 10.^(AmdB_rate_Lhigh/20);
    Y(m) = sqrt(sum(Am_rate_Lhigh.^2 .* h(m,1:Lhigh-1)));
    YdB(m) = 20*log10(Y(m));
  end
endfunction

% synth time domain waveform, broken into three frequency bands
function s = synth_time(Wo, L, Am, phase, N)
  s = zeros(1,N);
  t=0:N-1;
  for m=1:L
    s += Am(m)*exp(j*(Wo*m*t + phase(m)));
  end
endfunction

