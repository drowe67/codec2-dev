% plphase2.m
%
% Experiment to compare original and phase0 synthesised phase combined with rate K amplitude model.
% 1. Attempt to plot phase and group delay by removing linear phase component (n0), this doesn't
%    work very well for original phases, more meaningful for phase0.
% 2. Plot time domain speech in 0-1000Hz, 1000-2000, and 2000-4000Hz bands so we can explore envelope
%    of synthesised speech, which is realted to phase spectra in formants

#{
  Usage:

    $ cd codec2/build_linux
    $ ./src/c2sim ../raw/hts1a.raw --phase0 --dump hts1a
    $ ./src/c2sim ../raw/hts1a.raw --modelout - | ./misc/est_n0 > hts1a_n0.txt

    $ cd codec2/build_linux/octave
    $ octave-cli
    octave:> plphase2("../build_linux/hts1a", 44)
#}

function plphase2(samname, f, Nb=20, K=30)
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
  n0_name = strcat(samname,"_n0.txt");
  if (file_in_path(".",n0_name))
    n0 = load(n0_name);
  endif

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

  k = ' '; plot_group_delay=0; Pms = 6; plot_synth_sn=1; phase0_en=0;
  postfilter_en = 0;
  do
    Wo = model(f,1); F0 = Fs*Wo/(2*pi); L = model(f,2);
    Am = model(f,3:(L+2)); AmdB = 20*log10(Am);
    Am_freqs_kHz = (1:L)*Wo*4/pi;

    % resample from rate L to rate Lhigh (both linearly spaced)
    AmdB_rate_Lhigh = interp1([0 Am_freqs_kHz 4], [0 AmdB 0],
                              rate_Lhigh_sample_freqs_kHz, "spline", "extrap");

    % Filter at rate Lhigh, y = F(R(a)). Note we filter in linear energy domain,
    % and Lhigh are linearly spaced
    Y = zeros(1,Lhigh-1); YdB = zeros(1,Lhigh-1);
    for m=1:Lhigh-1
      Am_rate_Lhigh = 10.^(AmdB_rate_Lhigh/20);
      Y(m) = sqrt(sum(Am_rate_Lhigh.^2 .* h(m,1:Lhigh-1)));
      YdB(m) = 20*log10(Y(m));
    end

    if postfilter_en
      % straight line fit to YdB to estimate spectral slope SdB
      w = 2*pi*rate_Lhigh_sample_freqs_kHz*1000/Fs;
      st = round(200/F0high);
      en = round(3700/F0high);
      [m b] = linreg(w(st:en),YdB(st:en),en-st+1);
      SdB = w*m+b; S = 10.^(SdB/20);

      Y_energy1 = sum(Y .^ 2);

      % remove slope and expand dynamic range
      YdB -= SdB;
      YdB *= 2.0;
      YdB += SdB;

      % normalise energy
      Y = 10 .^ (YdB/20);
      Y_energy2 = sum(Y .^ 2);
      Y *= sqrt(Y_energy1/Y_energy2);
      YdB =20*log10(Y);
    end

    % resample from rate Lhigh to rate L (both linearly spaced)

    AmdB_ = interp1([0 rate_Lhigh_sample_freqs_kHz 4], [0 YdB 0], Am_freqs_kHz, "spline", "extrap");
    Am_ = 10 .^ (AmdB_/20);

    % remove linear component from original phase
    phase_linear = exp(j*(1:L)*Wo*n0(f));
    phase_rect = exp(j*phase(f,1:L));
    phase_centred_rect = phase_rect .* conj(phase_linear);
    phase_centred = angle(phase_centred_rect);
    phase_centred = unwrap(phase_centred);

    % Synthesised phase0 model using Hilbert Transform
    Nfft=512;
    sample_freqs_kHz = (Fs/1000)*[0:Nfft/2]/Nfft;  % fft frequency grid (nonneg freqs)
    Gdbfk = interp1([0 rate_Lhigh_sample_freqs_kHz 4], [0 YdB 0], sample_freqs_kHz, "spline", "extrap");
    if postfilter_en
      Gdbfk *= 1.5;
    end
    [phase_ht s] = mag_to_phase(Gdbfk, Nfft);
    for m=1:L
      b = round(m*Wo*Nfft/(2*pi));
      phase0(f,m) = phase_ht(b);
    end

    % TODO phase from Codec 2 3200 ?  We know that sounds better, essentially
    % a phase model from time domain LPC

    % plot time domain speech ---------------------------------------------
    figure(1); clf;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    y_off = -5000;
    if phase0_en, phase_ratek(f,1:L) = phase0(f,1:L); else phase_ratek(f,1:L) = phase(f,1:L); end
    if plot_synth_sn
      N=length(s);
      s1_lo = zeros(1,N); s2_lo = zeros(1,N);
      s1_mid = y_off*ones(1,N); s2_mid = y_off*ones(1,N);
      s1_hi = 2*y_off*ones(1,N); s2_hi = 2*y_off*ones(1,N);
      t=0:N-1; f0 = Wo*Fs2/pi; P = Fs/f0;
      for m=1:round(L/4)
        s1_lo += Am(m)*cos(Wo*m*t + phase(f,m));
        s2_lo += Am_(m)*cos(Wo*m*t + phase_ratek(f,m));
      end
      for m=round(L/4)+1:round(L/2)
        s1_mid += Am(m)*cos(Wo*m*t + phase(f,m));
        s2_mid += Am_(m)*cos(Wo*m*t + phase_ratek(f,m));
      end
      for m=round(L/2)+1:L
        s1_hi += Am(m)*cos(Wo*m*t + phase(f,m));
        s2_hi += Am_(m)*cos(Wo*m*t + phase_ratek(f,m));
      end
      maxy = max([s1_lo]); miny = min([s1_hi]);
      maxy = ceil(maxy/5000)*5000; miny = floor(miny/5000)*5000;
      subplot(211); hold on; plot(s1_lo,'g'); plot(s1_mid,'r'); plot(s1_hi,'b'); hold off;
      axis([1 length(s) miny maxy]); grid; title('orig Am & orig phase');
      subplot(212); hold on; plot(s2_lo,'g'); plot(s2_mid,'r'); plot(s2_hi,'b'); hold off;
      if phase0_en, phase_str = 'phase0'; else phase_str = 'orig phase'; end
      axis([1 length(s) miny maxy]); grid; title(sprintf("Filtered Am & %s",phase_str));
    end
    if (k == 'p')
    endif

    figure(2); clf;
    plot((1:L)*Wo*4000/pi, 20*log10(Am),"g+-;Am;");
    hold on;
    plot(rate_Lhigh_sample_freqs_kHz*1000, YdB, ';rate Lhigh YdB;b+-');

    % estimate group and phase delay and optionally plot ------------------------

    group_delay = [0 -((phase_centred(2:L) - phase_centred(1:L-1))/Wo)*1000/Fs];
    phase_delay = ( -phase_centred(1:L) ./ ((1:L)*Wo) )*1000/Fs;
    group_delay_phase0 = [0 -((phase0(f,2:L) - phase0(f,1:L-1))/Wo)*1000/Fs];
    phase_delay_phase0 = ( -phase0(f,1:L) ./ ((1:L)*Wo) )*1000/Fs;
    x_group = (0.5 + (1:L))*Wo*Fs2/pi;
    x_phase = (1:L)*Wo*Fs2/pi;
    if phase0_en
       if plot_group_delay == 1
          [ax h1 h2] = plotyy((0:255)*Fs2/256, Sw(f,:), x_group, group_delay_phase0);
       end
       if plot_group_delay == 2
          [ax h1 h2] = plotyy((0:255)*Fs2/256, Sw(f,:), x_phase, phase_delay_phase0);
       end
    else
       if plot_group_delay == 1
          [ax h1 h2] = plotyy((0:255)*Fs2/256, Sw(f,:), x_group, group_delay);
       end
      if plot_group_delay == 2
          [ax h1 h2] = plotyy((0:255)*Fs2/256, Sw(f,:), x_phase, phase_delay);
       end
    end
    if plot_group_delay
      axis(ax(1), [1 Fs2 -10 80]);
      axis(ax(2), [1 Fs2 -Pms Pms]);
      set(h2,'color','black');
      set(ax(2),'ycolor','black');
      ylabel(ax(1),'Amplitude (dB)');
      if plot_group_delay == 1
        ylabel(ax(2),'Group Delay (ms)');
      else
        ylabel(ax(2),'Phase Delay (ms)');
      end
    else
      plot((0:255)*Fs2/256, Sw(f,:));
      axis([1 Fs2 -10 80]);
      ylabel('Amplitude (dB)');
    end
    hold off; xlabel('Frequency (Hz)'); grid;

    % print tp EPS -------------------------------------------------------------
    if (k == 'p')
    endif

    figure(3); clf;
    subplot(211);
    adj = 0;
    if mean(phase_centred) < -2*pi
      adj = 2*pi;
    end
    adj_ratek = 0;
    if mean(phase0(f,1:L)) < -2*pi
      adj_ratek = 2*pi;
    end
    plot((1:L)*Wo*Fs2/pi, phase_centred+adj, "-og;phase;");
    hold on;
    plot((1:L)*Wo*Fs2/pi, phase0(f,1:L)+adj_ratek, "-or;phase0;");
    axis([0 Fs2 -2*pi 2*pi]);
    hold off;
    subplot(212);
    if plot_group_delay
      plot(x_group, group_delay, "-og;group delay;");
      hold on; plot(x_group, group_delay_phase0, "-or;group delay phase0;"); hold off;
      axis([1 Fs2 -Pms Pms]);
    else
      plot(x_phase, phase_delay, "-og;phase delay;");
      hold on; plot(x_group, phase_delay_phase0, "-or;phase delay phase0;"); hold off;
      axis([1 Fs2 -Pms Pms]);
    end

    if (k == 'p')
    endif

    % interactive menu

    if plot_group_delay==0, s1="group/phase dly"; end
    if plot_group_delay==1, s1="[group]/phase dly"; end
    if plot_group_delay==2, s1="group/[phase] dly"; end
    if phase0_en; s2="orig/[phase0]"; else s2="[orig]/phase0"; end
    printf("\rframe: %d  menu: n-next  b-back  g-%s 0-%s p-png f-postFilter[%d] q-quit ", f, s1, s2, postfilter_en);
    fflush(stdout);
    k = kbhit();
    if (k == 'n')
      f = f + 1;
    endif
    if (k == 'b')
      f = f - 1;
    endif
    if (k == 'g')
      if plot_group_delay, plot_group_delay = mod(plot_group_delay+1,3); end
    end
    if k == '0',
      if phase0_en, phase0_en = 0; else phase0_en = 1; end
    end
    if k == 'f',
      if postfilter_en, postfilter_en = 0; else postfilter_en = 1; end
    end
  until (k == 'q')
  printf("\n");

endfunction
