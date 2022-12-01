% snr_curves_plot.m
%
% Script to plot data from unittest/snr_cuvres

1;

function snr_scatter(mode, fmt)
  snrch = load(sprintf("../unittest/snrch_%s.txt",mode));
  snroffset = load(sprintf("../unittest/snroffset_%s.txt",mode));
  snrch -= snroffset;
  for i=1:length(snrch)
    fn = sprintf('../unittest/snrest_%s_%d.txt',mode,i);
    if exist(fn,'file') == 2
      snrest=load(fn);
      if i == 1
        plot(snrch(i)*ones(1,length(snrest)), snrest, sprintf('%s;%s;',fmt,mode));
      else
        plot(snrch(i)*ones(1,length(snrest)), snrest, fmt);
     end
    end
  end
endfunction

figure(1); clf; hold on;
snr_scatter('datac0','b+')
snr_scatter('datac1','g+')
snr_scatter('datac3','r+')
xlabel('SNR (dB)'); ylabel('SNRest (dB)')
%legend({'datac0','datac3'})
