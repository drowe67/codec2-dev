% ratek_resampler_plot.m
% Support for plotting curves from ratek_resampler.sh

function ratek_resampler_plot(plot_fn, varargin)
  more off;
  figure(1); clf;
  hold on;
  
  i = 1; bits_offset = 0;
  while i<=length(varargin)
    if strcmp(varargin{i},"continue")
      % used for 2nd/3rd stage, continue from last file loaded
      bits_offset = log2(max(x(:,1))); i++;
    else
      fn = varargin{i}; i++;
      leg = varargin{i}; i++;
      x = load(fn);
      plot(bits_offset+log2(x(:,1)),x(:,2),leg);
      bits_offset = 0;
    end  
  end

  % plot variance against bits for decorrelated scalars,
  % e.g. 6dB/bit/VQ element, or each bit reduces variance by factor 4
  K = 30; bits=1:20; m = 6/30; var0 = 20; % var0 arbitrary choice
  scalar = var0*10.^(-m*bits/10)
  plot(bits,scalar,'b+-;scalar;');
  
  xlabel('bits'); ylabel('var dB*dB'); grid('minor');
  axis([0 20 0 20]);
  print("-dpng", plot_fn);
endfunction
