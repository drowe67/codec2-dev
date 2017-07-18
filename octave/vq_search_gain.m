%----------------------------------------------------------------------
% abs() search with a linear gain term

function [idx contrib errors test_ g mg sl] = vq_search_gain(vq, data)
  
  [nVec nCols] = size(vq);
  nRows = rows(data);

  error = zeros(1,nVec);
  g = zeros(nRows, nVec);
  diff  = zeros(nVec, nCols);
  errors = zeros(1, nRows);
  idx = zeros(1, nRows);
  contrib = zeros(nRows, nCols);
  test_ = zeros(nVec, nCols);
  weights = ones(1,nCols);
  %weights(nCols-10+1:nCols) = 0.25;

  for f=1:nRows
    target = data(f,:);
    for i=1:nVec
      % work out gain for best match

      g(f, i) = (sum(target) - sum(vq(i,:)))/nCols;
      diff(i,:) = target - vq(i,:) - g(f, i);
      diff(i,:) .* weights;

      % abs in dB is MSE in linear

      error(i) = mean(abs(diff(i,:)));

      %printf("f: %d i: %d g: %f error: %f\n", f, i, g(f, i), error(i));
    end
    [mn min_ind] = min(error);

    idx(f) = min_ind;

    errors(f) = mn; 
    contrib(f,:) = test_(f,:) = vq(min_ind,:) + g(f,min_ind);
  end
  mg = 1; sl = 0;
endfunction

