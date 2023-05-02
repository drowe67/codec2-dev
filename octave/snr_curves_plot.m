% snr_curves_plot.m
%
% Companion script for unittest/raw_data_curves

1;

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

function [snr_ch per] = snr_scatter(source, mode, channel, colour)
  suffix = sprintf("_%s_%s_%s",source, mode, channel);
  snr = load(sprintf("snr%s.txt",suffix));
  offset = load(sprintf("offset%s.txt",suffix));
  snr -= offset;
  snr_x = []; snrest_y = [];
  for i=1:length(snr)
    fn = sprintf('snrest%s_%d.txt',suffix,i);
    if exist(fn,'file') == 2
      snrest=load(fn);
      if i == length(snr)
        plot(snr(i)*ones(1,length(snrest)), snrest, sprintf('%s;%s %s;',colour,source,mode));
      else
        plot(snr(i)*ones(1,length(snrest)), snrest, sprintf('%s',colour));
      end
      snr_x = [snr_x snr(i)]; snrest_y = [snrest_y mean(snrest)];
    end
  end
  plot(snr_x, snrest_y, sprintf('%s', colour));
endfunction

function [snr_ch per] = per_snr(mode, colour)
  snrch = load(sprintf("snrch_%s.txt",mode));
  snroffset = load(sprintf("snroffset_%s.txt",mode));
  snrch -= snroffset; 
  per = load(sprintf("per_%s.txt",mode));
  plot(snrch, per, sprintf('%so-;%s;', colour, mode));
endfunction

function snrest_snr_screen(source, channel)
  clf; hold on;
  snr_scatter(source, 'datac0', channel,'b+-')
  snr_scatter(source, 'datac1', channel,'g+-')
  snr_scatter(source, 'datac3', channel,'r+-')
  snr_scatter(source, 'datac4', channel,'c+-')
  snr_scatter(source, 'datac13', channel,'m+-')
  xlabel('SNR (dB)'); ylabel('SNRest (dB)'); grid('minor');
  axis([-12 12 -12 12]);
  a = axis;
  plot([a(1) a(2)],[a(1) a(2)],'bk-');
  hold off; grid;
  if strcmp(source,'ctx')
    title(sprintf('SNR estimate versus SNR (%s) (no compression)', channel));
  else
    title(sprintf('SNR estimate versus SNR (%s) (with compression)', channel));
  end
  legend('location','northwest');
endfunction

function snrest_snr_print(source, channel)
  state_vec = set_graphics_state_print();
  snrest_snr_screen(source, channel);
  print(sprintf("snrest_snr_%s.png", source), "-dpng", "-S1000,800");
  set_graphics_state_screen(state_vec);
endfunction

function ber_per_v_snr(source, mode, channel, colour)
  suffix = sprintf("_%s_%s_%s.txt",source, mode, channel);
  snr = load(sprintf("snr%s",suffix));
  offset = load(sprintf("offset%s",suffix));
  snr -= offset;
  ber = load(sprintf("ber%s",suffix)) + 1E-6;
  per = load(sprintf("per%s",suffix)) + 1E-6;
  semilogy(snr, ber, sprintf('%s;%s %s ber;', colour, source, mode));
  semilogy(snr, per, sprintf('%s;%s %s per;', colour, source, mode),'linewidth',3,'markersize',10);
endfunction

function per_v_snr(source, mode, channel, colour)
  suffix = sprintf("_%s_%s_%s.txt",source, mode, channel);
  snr = load(sprintf("snr%s",suffix));
  offset = load(sprintf("offset%s",suffix));
  snr -= offset;
  per = load(sprintf("per%s",suffix)) + 1E-6;
  if strcmp(channel,"awgn")
    semilogy(snr, per, sprintf('%s;%s %s;', colour, mode, channel));
  else
    semilogy(snr, per, sprintf('%s;%s %s;', colour, mode, channel),'linewidth',3,'markersize',10);
  end
endfunction

function thruput_v_snr(source, mode, channel, colour)
  suffix = sprintf("_%s_%s_%s.txt",source, mode, channel);
  snr = load(sprintf("snr%s",suffix));
  offset = load(sprintf("offset%s",suffix));
  snr -= offset;
  per = load(sprintf("per%s",suffix)) + 1E-6;
  if strcmp(mode,"datac0") Rb=291; end;
  if strcmp(mode,"datac1") Rb=980; end;
  if strcmp(mode,"datac3") Rb=321; end;
  if strcmp(mode,"datac4") Rb=87; end;
  if strcmp(mode,"datac13") Rb=65; end;
  if strcmp(channel,"awgn")
    plot(snr, Rb*(1-per), sprintf('%s;%s %s;', colour, mode, channel));
  else
    plot(snr, Rb*(1-per), sprintf('%s;%s %s;', colour, mode, channel),'linewidth',3,'markersize',10);
  end
endfunction

function octave_ch_noise_screen(channel)
  clf; hold on;
  ber_per_v_snr('oct','datac0',channel,'bo-')
  ber_per_v_snr('ch' ,'datac0',channel,'bx-')
  ber_per_v_snr('oct','datac1',channel,'go-')
  ber_per_v_snr('ch' ,'datac1',channel,'gx-')
  ber_per_v_snr('oct','datac3',channel,'ro-')
  ber_per_v_snr('ch' ,'datac3',channel,'rx-')
  xlabel('SNR (dB)'); grid;
  hold off;
  if strcmp(channel,"awgn")
    axis([-6 8 1E-3 1]);
  else
    axis([-2 12 1E-3 1]);
  end
  title(sprintf('Comparsion of Measuring SNR from Octave and ch tool (%s)', channel));
endfunction

function octave_ch_noise_print(channel)
  state_vec = set_graphics_state_print();
  octave_ch_noise_screen(channel);
  print(sprintf("octave_ch_noise_%s.png", channel), "-dpng","-S1000,800");
  set_graphics_state_screen(state_vec);
endfunction

function octave_c_tx_screen(channel)
  clf; hold on;
  ber_per_v_snr('oct','datac0',channel,'bo-')
  ber_per_v_snr('ctx','datac0',channel,'bx-')
  ber_per_v_snr('oct','datac1',channel,'go-')
  ber_per_v_snr('ctx','datac1',channel,'gx-')
  ber_per_v_snr('oct','datac3',channel,'ro-')
  ber_per_v_snr('ctx','datac3',channel,'rx-')
  xlabel('SNR (dB)'); grid;
  hold off;
  if strcmp(channel,"awgn")
    axis([-6 8 1E-3 1]);
  else
    axis([-2 12 1E-3 1]);
  end
  title(sprintf('Comparsion of Octave Tx and C Tx (no compression) (%s)', channel));
endfunction

function octave_c_tx_print(channel)
  state_vec = set_graphics_state_print();
  octave_c_tx_screen(channel);
  print(sprintf("octave_c_tx_%s.png", channel), "-dpng","-S1000,800");
  set_graphics_state_screen(state_vec);
endfunction

function octave_c_tx_comp_screen(channel)
  clf; hold on;
  ber_per_v_snr('oct','datac0',channel,'bo-')
  ber_per_v_snr('ctxc','datac0',channel,'bx-')
  ber_per_v_snr('oct','datac1',channel,'go-')
  ber_per_v_snr('ctxc','datac1',channel,'gx-')
  ber_per_v_snr('oct','datac3',channel,'ro-')
  ber_per_v_snr('ctxc','datac3',channel,'rx-')
  xlabel('SNR (dB)'); grid;
  hold off;
  if strcmp(channel,"awgn")
    axis([-6 8 1E-3 1]);
  else
    axis([-2 12 1E-3 1]);
  end
  title(sprintf('Comparsion of Octave Tx and C Tx (with compression) (%s)', channel));
endfunction

function octave_c_tx_comp_print(channel)
  state_vec = set_graphics_state_print();
  octave_c_tx_comp_screen(channel);
  print(sprintf("octave_c_tx_comp_%s.png", channel), "-dpng","-S1000,800");
  set_graphics_state_screen(state_vec);
endfunction

% composite AWGN and MPP for compressed
function c_tx_comp_screen
  clf; hold on;
  per_v_snr('ctxc','datac0','awgn','bo-')
  per_v_snr('ctxc','datac1','awgn','go-')
  per_v_snr('ctxc','datac3','awgn','ro-')
  per_v_snr('ctxc','datac4','awgn','co-')
  per_v_snr('ctxc','datac13','awgn','mo-')
  per_v_snr('ctxc','datac0','mpp','bx-')
  per_v_snr('ctxc','datac1','mpp','gx-')
  per_v_snr('ctxc','datac3','mpp','rx-')
  per_v_snr('ctxc','datac4','mpp','cx-')
  per_v_snr('ctxc','datac13','mpp','mx-')
  xlabel('SNR (dB)'); ylabel('PER'); grid;
  hold off;
  axis([-10 10 1E-3 1]);
  title('PER of C Raw Data Modes (with compression)');
endfunction

function c_tx_comp_print;
  state_vec = set_graphics_state_print();
  c_tx_comp_screen;
  print("c_tx_comp.png", "-dpng","-S1000,800");
  set_graphics_state_screen(state_vec);
endfunction

function c_tx_comp_thruput_screen
  clf; hold on;
  thruput_v_snr('ctxc','datac0','awgn','bo-')
  thruput_v_snr('ctxc','datac1','awgn','go-')
  thruput_v_snr('ctxc','datac3','awgn','ro-')
  thruput_v_snr('ctxc','datac4','awgn','co-')
  thruput_v_snr('ctxc','datac13','awgn','mo-')
  thruput_v_snr('ctxc','datac0','mpp','bx-')
  thruput_v_snr('ctxc','datac1','mpp','gx-')
  thruput_v_snr('ctxc','datac3','mpp','rx-')
  thruput_v_snr('ctxc','datac4','mpp','cx-')
  thruput_v_snr('ctxc','datac13','mpp','mx-')
  xlabel('SNR (dB)'); ylabel('bits/s'); grid;
  hold off;
  axis([-10 10 0 1000]);
  title(' Throughput for C Tx (with compression)');
  legend('location','west');
endfunction

function c_tx_comp_thruput_print;
  state_vec = set_graphics_state_print;
  c_tx_comp_thruput_screen;
  print("c_tx_comp_thruput.png", "-dpng","-S1000,800");
  set_graphics_state_screen(state_vec);
endfunction

#{
figure(1); octave_ch_noise_screen;
figure(2); octave_c_tx_screen;
figure(3); octave_c_tx_comp_screen
figure(4); snrest_snr_screen;

figure(5); octave_ch_noise_print;
figure(6); octave_c_tx_print;
figure(7); octave_c_tx_comp_print;
figure(8); snrest_snr_print;
#}
