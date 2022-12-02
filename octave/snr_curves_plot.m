% snr_curves_plot.m
%
% Script to plot data from unittest/snr_cuvres

1;

function [snr_ch per] = snr_scatter(mode, colour)
  snrch = load(sprintf("../unittest/snrch_%s.txt",mode));
  snroffset = load(sprintf("../unittest/snroffset_%s.txt",mode));
  snrch -= snroffset; snrch_x = []; snrest_y = [];
  for i=1:length(snrch)
    fn = sprintf('../unittest/snrest_%s_%d.txt',mode,i);
    if exist(fn,'file') == 2
      snrest=load(fn);
      if i == 1
        plot(snrch(i)*ones(1,length(snrest)), snrest, sprintf('%s+;%s;',colour,mode));
      else
        plot(snrch(i)*ones(1,length(snrest)), snrest, sprintf('%s+',colour));
      end
      snrch_x = [snrch_x snrch(i)]; snrest_y = [snrest_y mean(snrest)];
    end
  end
  plot(snrch_x, snrest_y, sprintf('%so-', colour));
endfunction

function [snr_ch per] = per_snr(mode, colour)
  snrch = load(sprintf("../unittest/snrch_%s.txt",mode));
  snroffset = load(sprintf("../unittest/snroffset_%s.txt",mode));
  snrch -= snroffset; 
  per = load(sprintf("../unittest/per_%s.txt",mode));
  plot(snrch, per, sprintf('%so-;%s;', colour, mode));
endfunction

% we need different font sizes for printing
function snr_scatter_screen
  clf; hold on;
  snr_scatter('datac0','b')
  snr_scatter('datac1','g')
  snr_scatter('datac3','r')
  xlabel('SNR (dB)'); ylabel('SNRest (dB)'); grid('minor');
  a = axis;
  plot([a(1) a(2)],[a(1) a(2)]);
  hold off;
endfunction

% we need different font sizes for printing
function snr_scatter_print(png_name)
  textfontsize = get(0,"defaulttextfontsize");
  linewidth = get(0,"defaultlinelinewidth");
  markersize = get(0, "defaultlinemarkersize");
  set(0, "defaulttextfontsize", 10);
  set(0, "defaultaxesfontsize", 10);
  set(0, "defaultlinelinewidth", 0.5);
  
  snr_scatter_screen;
  print("snr_curves.png", "-dpng", "-S500,500");

  % restore plot defaults
  set(0, "defaulttextfontsize", textfontsize);
  set(0, "defaultaxesfontsize", textfontsize);
  set(0, "defaultlinelinewidth", linewidth);  
  set(0, "defaultlinemarkersize", markersize);
endfunction

% we need different font sizes for printing
function per_snr_screen
  clf; hold on;
  per_snr('datac0','b')
  per_snr('datac1','g')
  per_snr('datac3','r')
  xlabel('SNR (dB)'); ylabel('PER'); grid;
  hold off;
endfunction

% we need different font sizes for printing
function per_snr_print
  textfontsize = get(0,"defaulttextfontsize");
  linewidth = get(0,"defaultlinelinewidth");
  markersize = get(0, "defaultlinemarkersize");
  set(0, "defaulttextfontsize", 10);
  set(0, "defaultaxesfontsize", 10);
  set(0, "defaultlinelinewidth", 0.5);
  
  per_snr_screen;
  print("per_snr.png", "-dpng", "-S500,500");

  % restore plot defaults
  set(0, "defaulttextfontsize", textfontsize);
  set(0, "defaultaxesfontsize", textfontsize);
  set(0, "defaultlinelinewidth", linewidth);  
  set(0, "defaultlinemarkersize", markersize);
endfunction

figure(1); snr_scatter_screen;
figure(2); snr_scatter_print;
figure(3); per_snr_screen;
figure(4); per_snr_print;

