%----------------------------------------------------------------------
% abs() search with a linear and slope term

function [idx contrib errors test_ g mg sl] = vq_search_slope(vq, data)
  [nVec nCols] = size(vq);
  nRows = rows(data);
  
  g = mg = sl = zeros(nRows, nVec);
  diff  = zeros(nVec, nCols);
  idx = errors = zeros(1, nRows);
  error = zeros(1, nVec);
  contrib = zeros(nRows, nCols);
  test_ = zeros(nVec, nCols);

  weights = ones(1,nCols);

  A  = [sum(1:nCols) (1:nCols)*(1:nCols)'; nCols  sum(1:nCols)];

  for f=1:nRows
    target = data(f,:);
    %target = 2*vq(1,:)+1;

    for i=1:nVec
      c = [target*(1:nCols)' sum(target)]';
      b = inv(A(:,:,i))*c;

      g(f,i) = b(1); sl(f,i) = b(2);
      
      diff(i,:) = target - (vq(i,:) + g(f,i) + sl(f,i)*(1:nCols));
      diff(i,:) .* weights;

      % abs in dB is MSE in linear

      error(i) = mean(abs(diff(i,:)));

      %printf("f: %d i: %d mg: %f g: %f sl: %f error: %f\n", f, i, mg(f,i), g(f,i), sl(f,i), error(i));
    end
    
    [mn min_ind] = min(error);
    errors(f) = mn; 
    idx(f) = min_ind(1);
   
    contrib(f,:) = test_(f,:) = vq(min_ind,:) + g(f,min_ind) + sl(f,min_ind)*(1:nCols);
  end

endfunction
