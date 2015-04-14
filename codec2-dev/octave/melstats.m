% melstats.m
% David Rowe April 2015
% 
% plots some stats of mel/lsp quantisers

function melstats(filename)

  mel = load(filename);
  [m n] = size(mel);
  nbins = 10;

  % histograms of each value

  figure(1)
  clf
  subplot(211)
  [h x] = hist(mel(:,1),nbins);
  plot(x,h,"1");
  hold on

  for i=2:n
    [h x] = hist(mel(:,i),nbins);
    colour = sprintf("%d",i);
    plot(x,h,colour);
  end
  hold off

  % histograms differences

  subplot(212)
  [h x] = hist(mel(:,1),nbins);
  plot(x,h,"1");
  hold on

  for i=2:n
    [h x] = hist(mel(:,i)-mel(:,i-1),nbins);
    colour = sprintf("%d",i);
    plot(x,h, colour);
  end
  hold off

  figure(2)
  plot(mel(:,1),mel(:,2),'r+')
  hold on;
  plot(mel(:,3),mel(:,4),'g+')
  plot(mel(:,5),mel(:,6),'b+')
  hold off;

endfunction
