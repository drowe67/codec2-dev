% plphase2.m
%
% Experiment to explore original and phase0 synthesised phase combined with rate K amplitude model.
% 1. Plot time domain speech in 0-1000Hz, 1000-2000, and 2000-4000Hz bands so we can explore envelope
%    of synthesised speech, which is related to phase spectra in formants
% 2. UI to cycle between orig/phase0 phase, and apply amplitude and phase post filter
% 3. Plot phase and group delay for phase0.  In earlier versions of this script we tried to plot
%    for this orginal phases but we couldn't get it to work.

#{
  Usage:

    $ cd codec2/build_linux
    $ ./src/c2sim ../raw/hts1a.raw --phase0 --dump hts1a

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

  k = ' '; plot_group_delay=0; Pms = 6; plot_synth_sn=1; phase0_en=1;
  postfilter_en = 0; ratek_en = 1;
  do
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
      YdB = amplitude_postfilter(rate_Lhigh_sample_freqs_kHz, YdB, Fs, F0high);
    end

    % Synthesised phase0 model using Hilbert Transform
    phase0(f,1:L) = synth_phase_from_mag(rate_Lhigh_sample_freqs_kHz, YdB, Fs, Wo, L, postfilter_en);

    % resample from rate Lhigh to rate L (both linearly spaced)
    AmdB_ = interp1([0 rate_Lhigh_sample_freqs_kHz 4], [0 YdB 0], Am_freqs_kHz, "spline", "extrap");
    Am_ = 10 .^ (AmdB_/20);

    % plot time domain speech ---------------------------------------------
    figure(1); clf;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    y_off = -5000;
    if phase0_en, phase_ratek(f,1:L) = phase0(f,1:L); else phase_ratek(f,1:L) = phase(f,1:L); end
    if plot_synth_sn
      N=length(s);
      [s1_lo s1_mid s1_hi] = synth_time(Wo, L, Am, phase(f,:), N);
      [s2_lo s2_mid s2_hi] = synth_time(Wo, L, Am_, phase_ratek(f,:), N);
      maxy =  10000;
      miny = -15000;
      maxy = ceil(maxy/5000)*5000; miny = floor(miny/5000)*5000;
      subplot(211); plot_time(s1_lo, s1_mid, s1_hi);
      axis([1 length(s) miny maxy]); grid; title('orig Am and orig phase');
      subplot(212); plot_time(s2_lo, s2_mid, s2_hi);
      if ratek_en, am_str = "filtered Am"; else am_str = "orig Am"; end
      if phase0_en, phase_str = 'phase0'; else phase_str = 'orig phase'; end
      if postfilter_en, phase_str = sprintf('%s and postfilter', phase_str); end
      axis([1 length(s) miny maxy]); grid; title(sprintf("%s and %s",am_str, phase_str));
    end

    figure(2); clf;
    if postfilter_en == 0
      plot((1:L)*Wo*4000/pi, 20*log10(Am),"g+-;Am;");
    endif
    hold on;
    if postfilter_en
      plot(rate_Lhigh_sample_freqs_kHz*1000, YdB_orig, ';rate Lhigh YdB;b-');
      plot(rate_Lhigh_sample_freqs_kHz*1000, YdB, ';rate Lhigh YdB postfilter;r-');
    else
      plot(rate_Lhigh_sample_freqs_kHz*1000, YdB, ';rate Lhigh YdB;b-');
    end

    % estimate group and phase delay and optionally plot ------------------------

    group_delay_phase0 = [0 -((phase0(f,2:L) - phase0(f,1:L-1))/Wo)*1000/Fs];
    phase_delay_phase0 = ( -phase0(f,1:L) ./ ((1:L)*Wo) )*1000/Fs;
    x_group = (0.5 + (1:L))*Wo*Fs2/pi;
    x_phase = (1:L)*Wo*Fs2/pi;
    if plot_group_delay
       if plot_group_delay == 1
          [ax h1 h2] = plotyy((0:255)*Fs2/256, Sw(f,:), x_group, group_delay_phase0);
       end
       if plot_group_delay == 2
          [ax h1 h2] = plotyy((0:255)*Fs2/256, Sw(f,:), x_phase, phase_delay_phase0);
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
      end
    else
      plot((0:255)*Fs2/256, Sw(f,:));
      axis([1 Fs2 -10 80]);
      ylabel('Amplitude (dB)');
      legend("boxoff"); legend("location","north");
    end
    hold off; xlabel('Frequency (Hz)');
    grid;

    figure(3); clf;
    subplot(211);
    adj_ratek = 0;
    if mean(phase0(f,1:L)) < -2*pi
      adj_ratek = 2*pi;
    end
    plot((1:L)*Wo*Fs2/pi, phase0(f,1:L)+adj_ratek, "-or;phase0;");
    axis([0 Fs2 -2*pi 2*pi]); ylabel('Phase (rads)')
    subplot(212);
    if plot_group_delay == 1
      plot(x_group, group_delay_phase0, "-or;group delay phase0;");
      axis([1 Fs2 -Pms Pms]); ylabel('Group delay (ms)')
    end
    if plot_group_delay == 2
      plot(x_group, phase_delay_phase0, "-or;phase delay phase0;");
      axis([1 Fs2 -Pms Pms]); ylabel('Phase delay (ms)')
    end

    if (k == 'p')
      [dir name ext]=fileparts(samname);

      % orig Am, orig phase
      [s1_lo s1_mid s1_hi] = synth_time(Wo, L, Am, phase(f,:), N);

      % filtered Am, phase0
      [YdB Y] = filter_rate_Lhigh(Lhigh,h,AmdB_rate_Lhigh);
      phase0 = synth_phase_from_mag(rate_Lhigh_sample_freqs_kHz, YdB, Fs, Wo, L, 0);
      AmdB_ = interp1([0 rate_Lhigh_sample_freqs_kHz 4], [0 YdB 0], Am_freqs_kHz, "spline", "extrap");
      Am_ = 10 .^ (AmdB_/20);
      [s2_lo s2_mid s2_hi] = synth_time(Wo, L, Am_, phase0, N);

      % filtered Am, phase0, postfilter
      YdB = amplitude_postfilter(rate_Lhigh_sample_freqs_kHz, YdB, Fs, F0high);
      phase0 = synth_phase_from_mag(rate_Lhigh_sample_freqs_kHz, YdB, Fs, Wo, L, 1);
      AmdB_ = interp1([0 rate_Lhigh_sample_freqs_kHz 4], [0 YdB 0], Am_freqs_kHz, "spline", "extrap");
      Am_ = 10 .^ (AmdB_/20);
      [s3_lo s3_mid s3_hi] = synth_time(Wo, L, Am_, phase0, N);

      figure(4); clf;
      maxy =  10000; miny = -15000;
      subplot(311); plot_time(s1_lo, s1_mid, s1_hi);
      axis([1 length(s) miny maxy]); grid; title('orig Am and orig phase');
      subplot(312); plot_time(s2_lo, s2_mid, s2_hi);
      axis([1 length(s) miny maxy]); grid; title('filtered Am and phase0');
      subplot(313); plot_time(s3_lo, s3_mid, s3_hi);
      axis([1 length(s) miny maxy]); grid; title('filtered Am and phase0 and post filter');

      fn=sprintf("plphase2_%s_%d_time",name,f);
      print(fn,"-depslatex","-S300,450");
      printf("\nprinting... %s\n", fn);

      figure(2);
      fn=sprintf("plphase2_%s_%d_freq",name,f);
      print(fn,"-depslatex","-S300,300");
      printf("printing... %s\n", fn);
    endif

    % interactive menu

    if plot_group_delay==0, s1="group/phase dly"; end
    if plot_group_delay==1, s1="[group]/phase dly"; end
    if plot_group_delay==2, s1="group/[phase] dly"; end
    if phase0_en; s2="orig/[phase0]"; else s2="[orig]/phase0"; end
    printf("\rframe: %d  n-nxt b-bk g-%s ratek-r 0-%s f-pf[%d] q-quit",f,s1,s2,postfilter_en);
    fflush(stdout);
    k = kbhit();
    if (k == 'n')
      f = f + 1;
    endif
    if (k == 'b')
      f = f - 1;
    endif
    if (k == 'g')
      plot_group_delay = mod(plot_group_delay+1,3);
    end
    if k == '0',
      if phase0_en, phase0_en = 0; else phase0_en = 1; end
    end
    if k == 'f',
      if postfilter_en, postfilter_en = 0; else postfilter_en = 1; end
    end
    if k == 'r',
      if ratek_en, ratek_en = 0; else ratek_en = 1; end
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
function [lo mid hi] = synth_time(Wo, L, Am, phase, N)
  lo = zeros(1,N);
  mid = ones(1,N);
  hi = ones(1,N);
  t=0:N-1;
  for m=1:round(L/4)
    lo += Am(m)*exp(j*(Wo*m*t + phase(m)));
  end
  for m=round(L/4)+1:round(L/2)
    mid += Am(m)*exp(j*(Wo*m*t + phase(m)));
  end
  for m=round(L/2)+1:L
    hi += Am(m)*exp(j*(Wo*m*t + phase(m)));
  end
endfunction

function plot_time(lo, mid, hi)
  y_off = -5000;
  hold on;
  plot(real(lo),sprintf('g;%s;',papr(lo)));
  plot(y_off + real(mid),sprintf('r;%s;',papr(mid)));
  plot(2*y_off + real(hi),sprintf('b;%s;',papr(hi)));
  legend("boxoff")
  hold off;
endfunction
