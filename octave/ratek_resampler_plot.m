% ratek_resampler_plot.m
% Support for plotting curves from ratek_resampler.sh

function ratek_resampler_plot(plot_fn, varargin)
  more off;
  epslatex_path = getenv("epslatex_path");
  if length(epslatex_path)
    textfontsize = get(0,"defaulttextfontsize");
    linewidth = get(0,"defaultlinelinewidth");
    set(0, "defaulttextfontsize", 10);
    set(0, "defaultaxesfontsize", 10);
    set(0, "defaultlinelinewidth", 0.5);
  end
  
  figure(1); clf;
  hold on;
  
  i = 1; bits_offset = 0; state='idle';
  while i<=length(varargin)
    if strcmp(state,'idle')
      if strcmp(varargin{i},"continue")
        % used for 2nd/3rd stage, continue from last file loaded
        bits_offset += log2(max(x(:,1)));
        i++;
        next_state = "cont";
      else
        fn = varargin{i}; i++;
        leg = varargin{i}; i++;
        x = load(fn);
        semilogy(log2(x(:,1)),x(:,2),leg);
        next_state = 'idle';
        bits_offset = 0;
      end
    end
   
    if strcmp(state,'cont')
      fn = varargin{i}; i++;
      leg = varargin{i}; i++;
      x = load(fn);
      semilogy(bits_offset+log2(x(:,1)),x(:,2),leg);
      next_state = 'idle';
    end
       
    state = next_state;
  end

  % plot variance against bits for decorrelated scalars,
  % e.g. 6dB/bit/VQ element, or each bit reduces variance by factor 4
  K = 20; bits=1:24; m = 6/30; var0 = 10; % var0 arbitrary choice
  scalar = var0*10.^(-m*bits/10);
  plot(bits,scalar,'b--;K=20 scalar;');
  plot([min(bits) max(bits)], [4 4],'b--;4 dB*dB;');
  
  xlabel('bits'); ylabel('Eq dB*dB'); grid('minor');
  axis([5 30 1 20]);
  if length(epslatex_path)
    legend("boxoff");
    old_dir=cd(epslatex_path);
    [dir name ext] = fileparts(plot_fn);
    fn = sprintf("%s.tex",name);
    print(fn, "-depslatex","-S400,300");
    printf("printing... %s%s\n", epslatex_path, fn);
    cd(old_dir);
    set(0, "defaulttextfontsize", textfontsize);
    set(0, "defaultaxesfontsize", textfontsize);
    set(0, "defaultlinelinewidth", linewidth);
  else
    print("-dpng", plot_fn);
  end  
endfunction
