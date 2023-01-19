% ratek_resampler_plot.m
% Support for plotting curves from ratek_resampler.sh

function ratek_resampler_plot(plot_fn, varargin)
  more off;
  figure(1); clf;
  hold on;
  
  i = 1;
  while i<=length(varargin)
    fn = varargin{i}; i++;
    leg = varargin{i}; i++;
    x = load(fn);
    plot(log2(x(:,1)),x(:,2),leg);
  end  
  xlabel('bits'); ylabel('var dB*dB'); grid;
  print("-dpng", plot_fn);
endfunction
