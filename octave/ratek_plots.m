% ratek1_plots.m
%
% David Rowe 2022
%
% Plots for Rate K report.
%


function ratek_plots()
  more off;

  newamp_700c;
  Fs = 8000;  K = 20;

  figure(1); clf; 
  f = 0:4000; plot(f, ftomel(f));
  grid;
  set(gca, 'FontSize', 16);
  xlabel('Frequency f (Hz)'); ylabel('mel(f)');
  print("ratek_mel_fhz","-dpng","-S500,500");

  figure(1); clf; 
  k = 1:K; plot(k, warp(k,K));
  grid;
  set(gca, 'FontSize', 16);
  xlabel('k'); ylabel('f (Hz)');
  print("warp_fhz_k","-dpng","-S500,500");

  figure(1); clf; 
  Fs=8000; f0=200; L=floor(Fs/(2*f0)); Nb=10;
  hold on;
  for m=1:L
    plot((1:L)*f0,generate_filter(m,f0,L,Nb),'o-');
  end
  hold off; grid;
  set(gca, 'FontSize', 16);
  xlabel('f (Hz)'); ylabel('h(m=f/F0)');
  print("filters_h","-dpng","-S500,500");

endfunction
