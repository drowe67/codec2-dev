% vq_pager.m
%
% Interactive Octave script to inspect vectors, e.g. for newamp1.m vector quantisation
%
% Usage:
%
%   octave:14> vq_pager(vq)

function vq_pager(vq)
  grid_sz = 4;

  [vq_rows vq_cols] = size(vq);
  
  r = 1;

  % Keyboard loop 

  k = ' ';
  do 
    figure(1); clf;
    n_plot = min(grid_sz*grid_sz,vq_rows-r+1)
    printf("r: %d n_plot: %d\n", r, n_plot);
    for i = 1:n_plot;
      subplot(grid_sz,grid_sz,i);
      plot(vq(r+i-1,:));
      axis([1 vq_cols -20 20]);
    end

    % interactive menu 

    printf("\rr: %d  menu: n-next  b-back  q-quit", r);
    fflush(stdout);
    k = kbhit();

    if k == 'n'
      r = r + grid_sz*grid_sz;
      r = min(r, vq_rows);
    endif
    if k == 'b'
      r = r - grid_sz*grid_sz;
      r = max(r,1);
    endif
  until (k == 'q')
  printf("\n");

endfunction

