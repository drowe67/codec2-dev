%----------------------------------------------------------------------
% abs() search with a linear, ampl scaling, and slope term

function [idx contrib errors b_log2] = vq_search_para(vq, data)
  [nVec nCols] = size(vq);
  nRows = rows(data);
  
  g = mg = sl = zeros(nRows, 1);
  diff  = zeros(nVec, nCols);
  idx = errors = zeros(1, nRows);
  error = zeros(1, nVec);
  contrib = zeros(nRows, nCols);
  test_ = zeros(nVec, nCols);

  weights = ones(1,nCols);

  A = zeros(4,4,nVec);
  K = nCols; k=(1:K); k2 = k.^2; 
  
  
  b_log = zeros(nVec, 4);
  b_log2 = [];
  
  for i=1:nVec
    v = vq(i,:);
    
    A(:,:,i) = [v*v'    k2*v'    k*v'    sum(v); ...
                k2*v'   k2*k2'   k*k2'   k*k';   ...                
                k*v'    k*k2'    k*k'    sum(k); ...
                sum(v)  k*k'     sum(k)  K       ];
  end
  
  for f=1:nRows
    t = data(f,:);

    for i=1:nVec
      v = vq(i,:);      
      c = [t*v' t*k2' t*k' sum(t)]';
      b = inv(A(:,:,i))*c;
      % b(1) = max(b(1),0.5);
      % b(2) = max(b(2),-0.2); b(2) = min(b(2),0.1);
      diff(i,:) = t - (b(1)*v + b(2)*k2 + b(3)*k + b(4));
      b_log(i,:) = b; 
      error(i) = diff(i,:) * diff(i,:)';

      %printf("f: %d i: %d e: %f b(1): %f b(2): %f b(3): %f b(4): %f\n", f, i, error(i), b(1), b(2), b(3), b(4));
    end
    
    [mn min_ind] = min(error);
    errors(f) = mn; 
    idx(f) = min_ind(1);
    b = b_log(min_ind,:);
    v = vq(min_ind,:);
    
    printf("f: %d i: %d b(1): %f b(2): %f b(3): %f b(4): %f\n", f, idx(f), b(1), b(2), b(3), b(4));
    contrib(f,:) = b(1)*v + b(2)*k2 + b(3)*k + b(4);
    b_log2 = [b_log2; b];
  end

endfunction
