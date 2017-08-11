%----------------------------------------------------------------------
% abs() search with a linear, ampl scaling, and slope term

function [idx contrib errors b_log2] = vq_search_slope(vq, data, closed_quant_fn, open_quant_fn)
  [nVec nCols] = size(vq);
  nRows = rows(data);
  
  diff  = zeros(nVec, nCols);
  idx = errors = zeros(1, nRows);
  error = zeros(1, nVec);
  contrib = zeros(nRows, nCols);
  b_log = zeros(nVec, 3);
  b_log2 = [];

  weights = ones(1,nCols);

  A = zeros(3,3,nVec);
  for i=1:nVec
    A(:,:,i) = [ sum(vq(i,:))       sum(1:nCols)         nCols;        ...
                 (1:nCols)*vq(i,:)' (1:nCols)*(1:nCols)' sum(1:nCols); ...
                 vq(i,:)*vq(i,:)'   (1:nCols)*vq(i,:)'   sum(vq(i,:))  ];
  end

  for f=1:nRows
    target = data(f,:);

    for i=1:nVec
      c = [sum(target) target*(1:nCols)' target*vq(i,:)' ]';
      b = inv(A(:,:,i))*c;
      if nargin >= 3;
        b = feval(closed_quant_fn,b);
      end
      b_log(i,:) = b; 
      
      diff(i,:) = target - (b(1)*vq(i,:) + b(2)*(1:nCols) + b(3));

      error(i) = diff(i,:) * diff(i,:)';

      %printf("f: %d i: %d mg: %f g: %f sl: %f error: %f\n", f, i, b(1), b(2), b(3), error(i));
    end
    
    [mn min_ind] = min(error);
    errors(f) = mn; 
    idx(f) = min_ind(1);
    b = b_log(min_ind,:);

    % optional quantisation performed after error min loop
    % - note we recalc gain to match energy after quantisation
    
    if nargin == 4
      b = feval(open_quant_fn, b);
      b(3) = (sum(target) - sum(b(1)*vq(min_ind,:) + b(2)*(1:nCols)))/nCols;    
    end

    printf("f: %d i: %d mg: %f sl: %f g: %f\n", f, idx(f), b(1), b(2), b(3));
    b_log2(f,:) = b;
    contrib(f,:) = b(1)*vq(min_ind,:) + b(2)*(1:nCols) + b(3);
  end

endfunction
