% standard mean squared error search

function [idx contrib errors test_ g mg sl] = vq_search_mse(vq, data)
  [nVec nCols] = size(vq);
  nRows = rows(data);
  
  error = zeros(1,nVec);
  errors = zeros(1, nRows);
  idx = zeros(1, nRows);
  contrib = zeros(nRows, nCols);
  test_ = zeros(nVec, nCols);

  for f=1:nRows
    target = data(f,:);
    for i=1:nVec
      diff = target - vq(i,:);
      error(i) = diff * diff';
    end
    [mn min_ind] = min(error);
    errors(f) = mn; idx(f) = min_ind; contrib(f,:) = vq(min_ind,:);
    test_(f,:) = vq(min_ind,:);
  end

  g = 0; mg = 1; sl = 0; % dummys for this function
endfunction

