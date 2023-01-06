% vq_pager.m
% Plot VQ entries in a grid to visualise codebook

function vq_pager(vq_fn, K=30)
  graphics_toolkit ("gnuplot");
  figure(1); clf;
  rows=4; cols=4;
  vq=load_f32(vq_fn,K);

  [n_left tmp] = size(vq);
  start_ind = 0;
  k = ' ';
  while n_left && k != 'q'
      n_plot = min(rows*cols, n_left);
      for i=1:n_plot
        c = mod(i,rows) + 1;
        r = floor(i/cols) + 1;
        subplot(rows,cols,i);
        plot(vq(start_ind + i,:));
      end
      n_left -= rows*cols;
      start_ind += rows*cols;
      printf("\rindex %d to %d q-quit", start_ind-rows*cols, start_ind-1);
      k = kbhit();
  end
end
