% plphase2.m
%
% Plot phase modelling information from dump files, compare two synthetic phase
% spectra, derived from two sets of {Am} magnitudes.

#{
  TODO - should we be running with origial phase here?
  Usage:

    $ cd codec2/build_linux
    $ ./src/c2sim ../raw/hts1a.raw --phase0 --dump hts1a
    $ ./src/c2sim ../raw/hts1a.raw --rateK --phase0 --dump hts1a_ratek

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

  k = ' '; plot_group_delay=1; Pms = 6; plot_orig=1; plot_synth_sn=1;
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

    % resample from rate Lhigh to rate L (both linearly spaced)

    AmdB_ = interp1([0 rate_Lhigh_sample_freqs_kHz 4], [0 YdB 0], Am_freqs_kHz, "spline", "extrap");
    Am_ = 10 .^ (AmdB_/20);

    % remove linear component from original phase
    phase_linear = exp(j*(1:L)*Wo*n0(f));
    phase_rect = exp(j*phase(f,1:L));
    phase_centred_rect = phase_rect .* conj(phase_linear);
    phase_centred = angle(phase_centred_rect);
    phase_centred = unwrap(phase_centred);

    % TODO phase from HT for synth speech, way to move thru orig phases,
    % ratek ampl + orig phase, orig amps + phase0
    Nfft=512;
    sample_freqs_kHz = (Fs/1000)*[0:Nfft/2]/Nfft;  % fft frequency grid (nonneg freqs)
    Gdbfk = interp1([0 rate_Lhigh_sample_freqs_kHz 4], [0 YdB 0], sample_freqs_kHz, "spline", "extrap");
    [phase_ht s] = mag_to_phase(Gdbfk, Nfft);
    for m=1:L
      b = round(m*Wo*Nfft/(2*pi));
      phase_ratek(f,m) = phase_ht(b);
    end

    % TODO phase from Codec 2 3200 ?  We know that sounds better

    % time domain speech ---------------------------------------------
    % TODO separate HF/LF (window/colour/toggle)
    figure(1); clf;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    y_off = -5000;
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
      maxy = max([s1_lo s2_lo]); miny = min([s1_hi s2_hi]);
      maxy = ceil(maxy/5000)*5000; miny = floor(miny/5000)*5000;
      subplot(211); hold on; plot(s1_lo,'g'); plot(s1_mid,'r'); plot(s1_hi,'b'); hold off;
      axis([1 length(s) miny maxy]); grid;
      subplot(212); hold on; plot(s2_lo,'g'); plot(s2_mid,'r'); plot(s2_hi,'b'); hold off;
      axis([1 length(s) miny maxy]); grid;
    end
    if (k == 'p')
    endif

    figure(2); clf;
    plot((1:L)*Wo*4000/pi, 20*log10(Am),"g+-;Am;");
    hold on;
    plot((1:L)*Wo*4000/pi, 20*log10(Am_),"r+-;Am*;");
    plot(rate_Lhigh_sample_freqs_kHz*1000, YdB, ';rate Lhigh YdB;b+-');

    % estimate group and phase delay

    group_delay = [0 -((phase_centred(2:L) - phase_centred(1:L-1))/Wo)*1000/Fs];
    phase_delay = ( -phase_centred(1:L) ./ ((1:L)*Wo) )*1000/Fs;
    group_delay_ratek = [0 -((phase_ratek(f,2:L) - phase_ratek(f,1:L-1))/Wo)*1000/Fs];
    phase_delay_ratek = ( -phase_ratek(f,1:L) ./ ((1:L)*Wo) )*1000/Fs;
    x_group = (0.5 + (1:L))*Wo*Fs2/pi;
    x_phase = (1:L)*Wo*Fs2/pi;
    if plot_orig
       if plot_group_delay
          [ax h1 h2] = plotyy((0:255)*Fs2/256, Sw(f,:), x_group, group_delay);
       else
          [ax h1 h2] = plotyy((0:255)*Fs2/256, Sw(f,:), x_phase, phase_delay);
       end
    else
       if plot_group_delay
          [ax h1 h2] = plotyy((0:255)*Fs2/256, Sw(f,:), x_group, group_delay_ratek);
       else
          [ax h1 h2] = plotyy((0:255)*Fs2/256, Sw(f,:), x_phase, phase_delay_ratek);
       end
    end
    hold off;
    axis(ax(1), [1 Fs2 -10 80]);
    axis(ax(2), [1 Fs2 -Pms Pms]);
    set(h2,'color','black');
    set(ax(2),'ycolor','black');
    xlabel('Frequency (Hz)');
    ylabel(ax(1),'Amplitude (dB)');
    if plot_group_delay
      ylabel(ax(2),'Group Delay (ms)');
    else
      ylabel(ax(2),'Phase Delay (ms)');
    end
    grid;

    if (k == 'p')
    endif

    figure(3); clf;
    subplot(211);
    adj = 0;
    if mean(phase_centred) < -2*pi
      adj = 2*pi;
    end
    adj_ratek = 0;
    if mean(phase_ratek(f,1:L)) < -2*pi
      adj_ratek = 2*pi;
    end
    plot((1:L)*Wo*Fs2/pi, phase_centred+adj, "-og;phase;");
    hold on;
    plot((1:L)*Wo*Fs2/pi, phase_ratek(f,1:L)+adj_ratek, "-or;phase ratek;");
    % axis([0 Fs2 -2*pi 2*pi]);
    hold off;
    subplot(212);
    if plot_group_delay
      plot(x_group, group_delay, "-og;group;");
      hold on; plot(x_group, group_delay_ratek, "-or;group ratek;"); hold off;
      % axis([1 Fs2 -1 Pms]);
    else
      plot(x_phase, phase_delay, "-og;phase;");
      hold on; plot(x_group, phase_delay_ratek, "-or;phase ratek;"); hold off;
      % axis([1 Fs2 -1 1]);
    end

    if (k == 'p')
    endif

    % interactive menu

    if plot_group_delay; s1="[group dly]/phase dly"; else s1="group dly/[phase dly]"; end
    if plot_orig; s2="[orig]/rateK"; else s2="orig/[rateK]"; end
    printf("\rframe: %d  menu: n-next  b-back  g-%s o-%s p-png  q-quit ", f, s1,s2);
    fflush(stdout);
    k = kbhit();
    if (k == 'n')
      f = f + 1;
    endif
    if (k == 'b')
      f = f - 1;
    endif
    if (k == 'g')
      if plot_group_delay
        plot_group_delay = 0;
      	Pms=1;
      else
        plot_group_delay = 1;
	      Pms=6;
      end
    endif
    if (k == 'o')
      if plot_orig
        plot_orig = 0;
      else
        plot_orig = 1;
      end
    endif

  until (k == 'q')
  printf("\n");

endfunction
