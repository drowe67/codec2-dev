% snr_curves_plot.m
%
% Script to plot data from unittest/snr_cuvres

1;

function snr_scatter(mode, fmt)
  snrch = load(sprintf("../unittest/snrch_%s.txt",mode));
  for i=1:length(snrch)
    fn = sprintf('../unittest/snrest_%s_%d.txt',mode,i);
    if exist(fn,'file') == 2
      snrest=load(fn);
      plot(snrch(i)*ones(1,length(snrest)), snrest,fmt);
    end
  end
endfunction

figure(1); clf; hold on;
snr_scatter('datac3','b+')
xlabel('SNR (dB)'); ylabel('SNRest (dB)')
legend({'datac3'})
