% autotest.m
% David Rowe Mar 2015
%
% Helper functions to plot output of C verson and difference between Octave and C versions

1;

function stem_sig_and_error(plotnum, subplotnum, sig, error, titlestr, axisvec)
  global no_plot_list;

  if find(no_plot_list == plotnum)
    return;
  end
  figure(plotnum)
  subplot(subplotnum)
  stem(sig,'g;Octave version;');
  hold on;
  stem(error,'r;Octave - C version (hopefully 0);');
  hold off;
  if nargin == 6
    axis(axisvec);
  end
  title(titlestr);
endfunction


function plot_sig_and_error(plotnum, subplotnum, sig, error, titlestr, axisvec)
  global no_plot_list;

  if find(no_plot_list == plotnum)
    return;
  end

  figure(plotnum)
  subplot(subplotnum)
  plot(sig,'g;Octave version;');
  hold on;
  plot(error,'r;Octave - C version (hopefully 0);');
  hold off;
  if nargin == 6
    axis(axisvec);
  end
  title(titlestr);
endfunction


function check(a, b, test_name)
  global passes;
  global fails;

  [m n] = size(a);
  printf("%s", test_name);
  for i=1:(25-length(test_name))
    printf(".");
  end
  printf(": ");  
  
  e = sum(abs(a - b))/n;
  if e < 1E-3
    printf("OK\n");
    passes++;
  else
    printf("FAIL\n");
    fails++;
  end
endfunction


