% rate7_high_fbf.m
%
% David Rowe 2022
%
% Rate K Experiment 7 - interactive Octave script to explore frame by frame
%                       operation of rate K resampling
%                     - Resampling rate Am to rate Lhigh and filter
%                     - this version plots original spectra
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README.md for
%    instructions.
%   ~/codec2-dev/build_linux$ ./src/c2sim ../raw/big_dog.raw --hpf --dump big_dog
%   $ cd ~/codec2-dev/octave
%   octave:14> ratek7_high_fbf("../build_linux/big_dog",60)


function ratek7_high_fbf(samname, f, Nb=20, K=30)
  more off;

  newamp_700c; melvq; pf_en = 0;
  Fs = 8000; resampler = 'spline'; Lhigh = 80;

  % load up text files dumped from c2sim ---------------------------------------

  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);
  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);
  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);
  rate_K_sample_freqs_kHz = mel_sample_freqs_kHz(K);

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

    Y = zeros(1,Lhigh-1); YdB = zeros(1,Lhigh-1);
    for m=1:Lhigh-1
      Am_rate_Lhigh = 10.^(AmdB_rate_Lhigh/20);
      Y(m) = sqrt(sum(Am_rate_Lhigh.^2 .* h(m,1:Lhigh-1)));
      YdB(m) = 20*log10(Y(m));
    end

    figure(3); clf;
    hold on;
    plot((0:255)*4000/256, Sw(f,:),";Sw;");
    l = sprintf(";rate %d AmdB;g+-", L);
    plot((1:L)*Wo*4000/pi, AmdB, l);

    % straight line fit to YdB to estimate spectral slope SdB
    w = 2*pi*rate_Lhigh_sample_freqs_kHz*1000/Fs;
    st = round(200/F0high);
    en = round(3700/F0high);
    [m b] = linreg(w(st:en),YdB(st:en),en-st+1);
    SdB = w*m+b; S = 10.^(SdB/20);
    plot(rate_Lhigh_sample_freqs_kHz*1000, SdB,';slope SdB;ro-');
    plot(rate_Lhigh_sample_freqs_kHz*1000, YdB-SdB,';flattened YdB;ro-');

    % fit LPC model to YS
    order = 10;
    [a G] = lpc_from_mag_spectrum(w, Y./S, order);
    gamma = 1.0;
    P = freqz(a.*(gamma.^(0:order)),1,w);
    PdB = 20*log10(abs(P));
    plot(rate_Lhigh_sample_freqs_kHz*1000, PdB,';Post Filter P;c-');

    if pf_en
      YdB(2:length(PdB)-1) = YdB(2:length(PdB)-1) + PdB(1:length(PdB)-2);
    end
    plot(rate_Lhigh_sample_freqs_kHz*1000, YdB, ';rate Lhigh YdB;b+-');

    axis([0 Fs/2 -10 80]);
    hold off;

    % interactive menu ------------------------------------------

    printf("\rframe: %d  menu: n-next  b-back  q-quit p-print", f);
    fflush(stdout);
    k = kbhit();

    if k == 'n'
      f = f + 1;
    endif
    if k == 'b'
      f = f - 1;
    endif
    if k == 'f'
      if pf_en, pf_en = 0; else pf_en = 1; end
    end
    if (k == 'p')
      [dir name ext]=fileparts(samname);
      hl = legend({"Spectrum","Rate AmdB" "Rate Lhigh YdB"}, "location", "northeast");
      legend("boxoff")
      xlabel('Freq (Hz)'); ylabel('Amplitude (dB)');
      print(sprintf("ratek7_%s_%d",name,f),"-depsc","-S300,300");
      printf("printing...")
    endif

  until (k == 'q')
  printf("\n");

endfunction
