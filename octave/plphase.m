% Copyright David Rowe 2009
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%
% Plot phase modelling information from dump files.

#{  
  Usage:

  1/ Contrived speech-like test signal:
  
  $ cd codec2/build_linux
  $ ./misc/timpulse --f0 100 --n0 1 --filter | ./src/c2sim - --rateK --modelout - | ./misc/est_n0 > imp_n0.txt
  $ ./misc/timpulse --f0 100 --n0 1 --filter | ./src/c2sim - --rateK --phase0 --dump imp
  octave:> plphase("../build_linux/imp", 20)

  2/ A real speech file:
  
  Generate n0 estimates for original and phase0:

    $ cd codec2/build_linux
    $ ./src/c2sim ../raw/hts1a.raw --rateK --modelout - | ./misc/est_n0 > hts1a_n0.txt
    $ ./src/c2sim ../raw/hts1a.raw --rateK --phase0 --modelout - | ./misc/est_n0 > hts1a_phase0_n0.txt

  Generate "dump files":
  
    $ ./src/c2sim ../raw/hts1a.raw --rateK --phase0 --dump hts1a
    
  octave:> plphase("../build_linux/hts1a", 44)
#}

function plphase(samname, f)

  Fs = 8000; Fs2 = Fs/2;
  
  sn_name = strcat(samname,"_sn.txt");
  Sn = load(sn_name);

  sw_name = strcat(samname,"_sw.txt");
  Sw = load(sw_name);

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);

  n0_name = strcat(samname,"_n0.txt");
  if (file_in_path(".",n0_name))
    n0 = load(n0_name);
  endif
  phase0_n0_name = strcat(samname,"_phase0_n0.txt");
  if (file_in_path(".",phase0_n0_name))
    phase0_n0 = load(phase0_n0_name);
  endif
  
  phase_name = strcat(samname,"_phase.txt");
  if (file_in_path(".",phase_name))
    phase = load(phase_name);
  endif

  phase_name_ = strcat(samname,"_phase_.txt");
  if (file_in_path(".",phase_name_))
    phase_ = load(phase_name_);
  endif

  k = ' '; plot_group_delay=0; Pms = 1; plot_orig=1;
  do
    figure(1); clf;
    s = [ Sn(2*f-1,:) Sn(2*f,:) ];
    plot(s);
    grid;
    axis([1 length(s) -20000 20000]);
    if (k == 'p')
       pngname = sprintf("%s_%d_sn",samname,f);
       png(pngname);
    endif

    figure(2); clf;
    Wo = model(f,1);
    L = model(f,2);
    Am = model(f,3:(L+2));
    plot((1:L)*Wo*4000/pi, 20*log10(Am),"g+-;Am;");
    hold on;
 
    % estimate group and phase delay

    if plot_orig
      phase_f = phase(f,1:L);
      phase_linear = exp(j*(1:L)*Wo*n0(f));
    else
      phase_f = phase_(f,1:L);
      phase_linear = exp(j*(1:L)*Wo*phase0_n0(f));
    end
    phase_rect = exp(j*phase_f);
    phase_centred_rect = phase_rect .* conj(phase_linear);
    phase_centred = angle(phase_centred_rect);
    group_delay = [0 -(angle(phase_centred_rect(2:L).*conj(phase_centred_rect(1:L-1)))/Wo)*1000/Fs];
    phase_delay = ( -phase_centred ./ ((1:L)*Wo) )*1000/Fs;
    x_group = (0.5 + (1:L))*Wo*Fs2/pi;
    x_phase = (1:L)*Wo*Fs2/pi;
    if plot_group_delay
      ax = plotyy((0:255)*Fs2/256, Sw(f,:), x_group, group_delay);
    else
      ax = plotyy((0:255)*Fs2/256, Sw(f,:), x_phase, phase_delay);
    end
    hold off;
    axis(ax(1), [1 Fs2 -10 80]);
    axis(ax(2), [1 Fs2 -Pms Pms]);
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
    plot((1:L)*Wo*Fs2/pi, phase_f, "-o;phase;");
    hold on;
    plot((1:L)*Wo*Fs2/pi, phase_centred, "-og;phase centered;"); axis([0 Fs2 -pi pi]);
    hold off;
    subplot(212);
    if plot_group_delay
      plot(x_group, group_delay, "-o;group;");
      axis([1 Fs2 -Pms Pms]);
    else
      plot(x_phase, phase_delay, "-o;phase;");
      axis([1 Fs2 -1 1]);
    end
    
    if (k == 'p')
      pngname = sprintf("%s_%d_phase",samname,f);
      png(pngname);
    endif
 
    % interactive menu

    if plot_group_delay; s1="[group dly]/phase dly"; else s1="group dly/[phase dly]"; end
    if plot_orig; s2="[orig]/phase0"; else s2="orig/[phase0]"; end
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
