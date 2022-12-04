% snr_curves_plot.m
%
% Companion script for unittest/Script to plot data from unittest/snr_curves.sh

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

function state_vec = set_graphics_state_print()
  textfontsize = get(0,"defaulttextfontsize");
  linewidth = get(0,"defaultlinelinewidth");
  markersize = get(0, "defaultlinemarkersize");
  set(0, "defaulttextfontsize", 16);
  set(0, "defaultaxesfontsize", 16);
  set(0, "defaultlinelinewidth", 1);
  state_vec = [textfontsize linewidth markersize];
endfunction

function set_graphics_state_screen(state_vec) 
  textfontsize = state_vec(1);
  linewidth = state_vec(2);
  markersize = state_vec(3);
  set(0, "defaulttextfontsize", textfontsize);
  set(0, "defaultaxesfontsize", textfontsize);
  set(0, "defaultlinelinewidth", linewidth);  
  set(0, "defaultlinemarkersize", markersize);
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

function ber_per_v_snr(source, mode, colour)
  suffix = sprintf("_%s_%s.txt",source, mode);
  snr = load(sprintf("../unittest/snr%s",suffix));
  offset = load(sprintf("../unittest/offset%s",suffix));
  snr -= offset;
  ber = load(sprintf("../unittest/ber%s",suffix)) + 1E-6;
  per = load(sprintf("../unittest/per%s",suffix)) + 1E-6;
  semilogy(snr, ber, sprintf('%s;%s %s ber;', colour, source, mode));
  semilogy(snr, per, sprintf('%s;%s %s per;', colour, source, mode));
 endfunction

function ber_per_v_snr_screen
  clf; hold on;
  ber_per_v_snr('oct','datac0','bo-')
  ber_per_v_snr('ch','datac0','bx-')
  ber_per_v_snr('oct','datac1','go-')
  ber_per_v_snr('ch','datac1','gx-')
  ber_per_v_snr('oct','datac3','ro-')
  ber_per_v_snr('ch','datac3','rx-') #}
  xlabel('SNR (dB)'); grid;
  hold off; axis([-6 8 1E-3 1]);
endfunction

function ber_per_v_snr_print
  state_vec = set_graphics_state_print()
  ber_per_v_snr_screen;
  print("ber_per_v_snr.png", "-dpng","-S800,600");
  set_graphics_state_screen(state_vec);
endfunction

figure(1); ber_per_v_snr_screen;
figure(2); ber_per_v_snr_print;

#{
figure(1); snr_scatter_screen;
figure(2); snr_scatter_print;
figure(3); per_snr_screen;
figure(4); per_snr_print;
#}

