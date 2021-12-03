% plphase2.m
%
% Plot phase modelling information from dump files, compare two synthetic phase 
% spectra, derived from two sets of {Am} magnitudes.

#{  
  Usage:

    $ cd codec2/build_linux
    $ ./src/c2sim ../raw/hts1a.raw --phase0 --dump hts1a
    $ ./src/c2sim ../raw/hts1a.raw --rateK --phase0 --dump hts1a_ratek

    octave:> plphase2("../build_linux/hts1a", 44)
#}

function plphase2(samname, f)

  Fs = 8000; Fs2 = Fs/2;
  
  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);

  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);

  qmodel_name = strcat(samname,"_ratek_qmodel.txt");
  qmodel = load(qmodel_name);

  % H is the synthesised phase without linear component
  H_name = strcat(samname,"_H.txt");
  if (file_in_path(".",H_name))
    phase = unwrap(load(H_name),pi,2);
  endif
  ratek_H_name = strcat(samname,"_ratek_H.txt");
  if (file_in_path(".",ratek_H_name))
    phase_ratek = unwrap(load(ratek_H_name),pi,2);
  endif
  
  k = ' '; plot_group_delay=1; Pms = 6; plot_orig=1; plot_synth_sn=1;
  do
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    Am_ = qmodel(f,3:(L+2));
    
    figure(1); clf;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    if plot_synth_sn
      N=length(s);
      s1 = zeros(1,N); s2 = zeros(1,N); t=0:N-1; f0 = Wo*Fs2/pi; P = Fs/f0;
      for m=1:L
        s1 += Am(m)*cos(Wo*m*t + phase(f,m));
        s2 += Am_(m)*cos(Wo*m*t + phase_ratek(f,m));
      end
      plot(s1,'g'); hold on; plot(s2,'r'); hold off;
    end
    grid;
    axis([1 length(s) -20000 20000]);
    if (k == 'p')
       pngname = sprintf("%s_%d_sn",samname,f);
       png(pngname);
    endif

    figure(2); clf;
    plot((1:L)*Wo*4000/pi, 20*log10(Am),"g+-;Am;");
    hold on;
    plot((1:L)*Wo*4000/pi, 20*log10(Am_),"r+-;Am*;");
 
    % estimate group and phase delay

    group_delay = [0 -((phase(f,2:L) - phase(f,1:L-1))/Wo)*1000/Fs];
    phase_delay = ( -phase(f,1:L) ./ ((1:L)*Wo) )*1000/Fs;
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
    xlabel('Frequency (Hz)');
    ylabel(ax(1),'Amplitude (dB)');
    if plot_group_delay
      ylabel(ax(2),'Group Delay (ms)');
    else
      ylabel(ax(2),'Phase Delay (ms)');
    end
    grid;
    
    if (k == 'p')
       pngname = sprintf("%s_%d_sw",samname,f);
       png(pngname);
    endif

    figure(3); clf;
    subplot(211);
    plot((1:L)*Wo*Fs2/pi, phase(f,1:L), "-og;phase;");
    hold on;
    plot((1:L)*Wo*Fs2/pi, phase_ratek(f,1:L), "-or;phase ratek;");
    axis([0 Fs2 -2*pi 2*pi]);
    hold off;
    subplot(212);
    if plot_group_delay
      plot(x_group, group_delay, "-og;group;");
      hold on; plot(x_group, group_delay_ratek, "-or;group ratek;"); hold off;
      axis([1 Fs2 -Pms Pms]);
    else
      plot(x_phase, phase_delay, "-og;phase;");
      hold on; plot(x_group, phase_delay_ratek, "-or;phase ratek;"); hold off;
      axis([1 Fs2 -1 1]);
    end
    
    if (k == 'p')
      pngname = sprintf("%s_%d_phase",samname,f);
      png(pngname);
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

    % optional print to PNG
    if (k == 'p')
       pngname = sprintf("%s_%d",samname,f);
       png(pngname);
    endif

  until (k == 'q')
  printf("\n");

endfunction
