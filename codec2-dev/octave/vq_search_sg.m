%----------------------------------------------------------------------
% abs() search with a linear and slope term

function [idx contrib errors b_log2] = vq_search_sg(vq, data)
  [nVec nCols] = size(vq);
  nRows = rows(data);
  
  idx = errors = zeros(1, nRows);
  error = zeros(1, nVec);
  contrib = zeros(nRows, nCols);

  b_log = zeros(nVec, 2);
  b_log2 = [];

  k = 1:nCols;
  
  for f=1:nRows
    t = data(f,:);

    for i=1:nVec
      v = vq(i,:);
      A  = [k*k' sum(k); sum(k) nCols];
      c = [(t*k'-v*k') (sum(t)-sum(v))]';
      b = inv(A)*c;
      b(1) = quantise([-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0], b(1));
      b_log(i,:) = b; 
      diff = t - (v + b(1)*k + b(2));
      error(i) = diff*diff';
      %printf("  i: %d error %f\n", i, error(i));
    end
    
    [mn min_ind] = min(error);
    errors(f) = mn; 
    idx(f) = min_ind(1);
    b = b_log(min_ind,:);
    v = vq(min_ind,:);
    printf("f: %d idx: %d b(1): %f b(2): %f\n", f, idx(f), b(1), b(2));
    contrib(f,:) = v + b(1)*k + b(2);
    b_log2 = [b_log2; b];
  end

endfunction
