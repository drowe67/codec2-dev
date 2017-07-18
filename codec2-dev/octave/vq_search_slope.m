%----------------------------------------------------------------------
% abs() search with a linear, ampl scaling, and slope term

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

  A = zeros(3,3,nVec);
  for i=1:nVec
    A(:,:,i) = [sum(vq(i,:))      nCols        sum(1:nCols); ...
                vq(i,:)*vq(i,:)'   sum(vq(i,:)) (1:nCols)*vq(i,:)'; ...
                (1:nCols)*vq(i,:)' sum(1:nCols) (1:nCols)*(1:nCols)'];
  end

  for f=1:nRows
    target = data(f,:);
    %target = 2*vq(1,:)+1;

    for i=1:nVec
      c = [sum(target) target*vq(i,:)' target*(1:nCols)']';
      b = inv(A(:,:,i))*c;

      %b(1) = quantise([0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0], b(1));
      %b(3) = quantise([-1.0 -0.5 0.0 0.5 1.0], b(3));
      %b(1) = quantise([0.8 1.2], b(1));
      %b(3) = quantise([-0.2 0.2], b(3));
      mg(f,i) = b(1); g(f,i) = b(2); sl(f,i) = b(3);
      
      diff(i,:) = target - (mg(f,i)*vq(i,:) + g(f,i) + sl(f,i)*(1:nCols));
      diff(i,:) .* weights;

      % abs in dB is MSE in linear

      error(i) = mean(abs(diff(i,:)));

      %printf("f: %d i: %d mg: %f g: %f sl: %f error: %f\n", f, i, mg(f,i), g(f,i), sl(f,i), error(i));
    end
    
    [mn min_ind] = min(error);
    errors(f) = mn; 
    idx(f) = min_ind(1);
   
    contrib(f,:) = test_(f,:) = mg(f,min_ind)*vq(min_ind,:) + g(f,min_ind) + sl(f,min_ind)*(1:nCols);
  end

endfunction
