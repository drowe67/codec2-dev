% kmeans_tests.m
%
% David Rowe June 2017
%
%
% Trying a few variations on the kmeans algorithm for quantisation of
% spectral envelope

1;

%----------------------------------------------------------------------
%
% User defined search functions
%
%----------------------------------------------------------------------

% standard mean squared error search

function [idx errors test_ g mg] = vq_search_mse(vq, data)
  [nVec nCols] = size(vq);
  nRows = length(data);
  
  error = zeros(1,nVec);
  errors = zeros(1, nRows);
  idx = zeros(1, nRows);
  test_ = zeros(nVec, nCols);

  for f=1:nRows
    target = data(f,:);
    for i=1:nVec
      diff = target - vq(i,:);
      error(i) = diff * diff';
    end
    [mn min_ind] = min(error);
    errors(f) = mn; idx(f) = min_ind;
    test_(f,:) = vq(min_ind,:);
  end

  g = mg = 1; % dummys for this function
endfunction


% abs() search with a linear gain term

function [idx errors test_ g mg] = vq_search_gain(vq, data)
  
  [nVec nCols] = size(vq);
  nRows = length(data);

  error = zeros(1,nVec);
  g = zeros(nRows, nVec);
  diff  = zeros(nVec, nCols);
  errors = zeros(1, nRows);
  idx = zeros(1, nRows);
  test_ = zeros(nVec, nCols);

  for f=1:nRows
    target = data(f,:);
    for i=1:nVec
      % work out gain for best match

      g(f, i) = (sum(target) - sum(vq(i,:)))/nCols;
      diff(i,:) = target - vq(i,:) - g(f, i);

      % abs in dB is MSE in linear

      error(i) = mean(abs(diff(i,:)));

      %printf("f: %d i: %d g: %f error: %f\n", f, i, g(f, i), error(i));
    end
    [mn min_ind] = min(error);
    errors(f) = mn; idx(f) = min_ind;
    test_(f,:) = vq(min_ind,:) + g(f,min_ind);
  end
  mg = 1;
endfunction


% abs() search with a linear plus ampl scaling term

function [idx errors test_ g mg] = vq_search_mag(vq, data)
  [nVec nCols] = size(vq);
  nRows = length(data);
  
  g = mg = zeros(nRows, nVec);
  diff  = zeros(nVec, nCols);
  errors = zeros(1, nRows);
  idx = error = zeros(1, nVec);
  mg_log = [];
  test_ = zeros(nVec, nCols);

  weights = ones(1,nCols);
  %weights(1) = 0;

  for f=1:nRows
    target = data(f,:);
    %target = 2*vq(1,:)+1;

    for i=1:nVec
      % work out gain and amp scaling for best match

      A = [sum(vq(i,:)) nCols; vq(i,:)*vq(i,:)' sum(vq(i,:))];
      c = [sum(target) target*vq(i,:)']';
      b = inv(A)*c;

      g(f,i) = b(2);  mg(f,i) = b(1);
      diff(i,:) = target - (mg(f,i)*vq(i,:) + g(f,i));
      diff(i,:) .* weights;

      % abs in dB is MSE in linear

      error(i) = mean(abs(diff(i,:)));

      %printf("f: %d i: %d mg: %f g: %f error: %f\n", f, i, mg(f,i), g(f,i), error(i));
    end
    
    [mn min_ind] = min(error);
    errors(f) = mn; 
    idx(f) = min_ind;
    test_(f,:) = mg(f,min_ind) * vq(min_ind,:) + g(f,min_ind);
    mg_log = [mg_log mg(f,min_ind)];
  end

  %figure(2); clf; hist(mg_log);
endfunction


% evaluate database test using vq, with selectable search function.  Can be operated in
% GUI mode to analyse in fine detail or batch mode to evaluate lots of data.

function sd_per_frame = run_test(vq, test, nVec, search_func, gui_en = 1)

  % Test VQ using test data -----------------------

  % Each test must output:
  %   errors: - error for best vector
  %   idx...: - index of best vector
  %   test_.: - best approximation of target test database
  
  [nRows nCols] = size(test);

  [idx errors test_ g mg] = feval(search_func, vq, test);
  
  % sd over time
 
  sd_per_frame = zeros(nRows,1);
  for i=1:nRows
    sd_per_frame(i) = std(test(i,:) - test_(i,:));
  end

  printf("  mode: %4s sd: %3.2f dB\n", search_func, mean(sd_per_frame));
  
  % plots sd and errors over time

  if gui_en
    figure(1); clf; subplot(211); plot(sd_per_frame); title('SD'); subplot(212); plot(errors); title('mean error');

    % display m frames, printing some stats, plotting vector to give visual idea of match

    figure(2); clf;
    [errors_dec frame_dec] = sort(errors, "descend");
    m = 4;
    for i=1:m
      af = frame_dec(i); aind = idx(af);
      l = sprintf("idx: %d", aind);
      ag = 0; amg = 1;
      if strcmp(search_func, "vq_search_gain") || strcmp(search_func, "vq_search_mag")
        ag = g(af,aind);
        l = sprintf("%s g: %3.2f", l, ag);
      end
      if strcmp(search_func, "vq_search_mag")
        amg = mg(af,aind);
        l = sprintf("%s mg: %3.2f", l, amg);
      end
      %printf("%d f: %d %s\n", i, af, l);

      subplot(sqrt(m),sqrt(m),i);
      l1 = sprintf("b-;fr %d;", af);
      plot(test(af,:), l1); 

      hold on; 
      l2 = sprintf("g-+;ind %d;", aind); 
      plot(vq(aind, :), l2); 
      l3 = sprintf("g-o;%s;",l);   
      plot(amg*vq(aind, :) + ag, l3); 
      hold off;
      axis([1 nCols -10 40]);
    end
  end  
  
endfunction


function compare_hist(atitle, sdpf_mse, sdpf_gain, sdpf_mag)
  [mse_yy, mse_xx] = hist(sdpf_mse);
  [gain_yy, gain_xx] = hist(sdpf_gain);
  [mag_yy, mag_xx] = hist(sdpf_mag);

  plot(mse_xx, mse_yy, 'b+-;mse;');
  hold on;
  plot(gain_xx, gain_yy, 'g+-;gain;');
  plot(mag_xx, mag_yy, 'r+-;mag;');
  hold off;
  title(atitle)
end


function long_tests(quick_check=0)

  K = 10;
  load surf_train_120;
  load surf_all;
  if quick_check
    NtrainVec = 1000;
    NtestVec = 100;
  else
    NtrainVec = length(surf_train_120);
    NtestVec = length(surf_all);
  end
  trainvec = surf_train_120(1:NtrainVec,1:K);
  testvec = surf_all(1:NtestVec,1:K);

  % Test 1 -------------------------------------------------------

  % standard kmeans, conventional MSE based training

  printf("Nvec = 64\n");

  Nvec = 64; [idx vq]  = kmeans(trainvec, Nvec, "emptyaction", "singleton");

  sdpf_mse  = run_test(vq, testvec, Nvec, 'vq_search_mse', gui_en=0);
  sdpf_gain = run_test(vq, testvec, Nvec, 'vq_search_gain', gui_en=0);
  sdpf_mag = run_test(vq, testvec, Nvec, 'vq_search_mag', gui_en=0);

  figure(1); clf; compare_hist("K=64 SD Histograms", sdpf_mse, sdpf_gain, sdpf_mag)

  % Test 2 --------------------------------------

  printf("Nvec = 256\n");

  Nvec = 256; [idx vq]  = kmeans(trainvec, Nvec, "emptyaction", "singleton");

  sdpf_mse  = run_test(vq, testvec, Nvec, 'vq_search_mse', gui_en=0);
  sdpf_gain = run_test(vq, testvec, Nvec, 'vq_search_gain', gui_en=0);
  sdpf_mag = run_test(vq, testvec, Nvec, 'vq_search_mag', gui_en=0);

  figure(2); clf; compare_hist("K=256 SD Histograms", sdpf_mse, sdpf_gain, sdpf_mag)

  % Test 3 --------------------------------------

  printf("Nvec = 64 150Hz HPF on train and test\n");

  load surf_train_120_hpf150;
  trainvec = surf_train_120_hpf150(1:NtrainVec,1:K);
  load surf_all_hpf150;
  testvec = surf_all_hpf150(1:NtestVec,1:K);

  Nvec = 64; [idx vq]  = kmeans(trainvec, Nvec, "emptyaction", "singleton");

  sdpf_mse  = run_test(vq, testvec, Nvec, 'vq_search_mse', gui_en=0);
  sdpf_gain = run_test(vq, testvec, Nvec, 'vq_search_gain', gui_en=0);
  sdpf_mag = run_test(vq, testvec, Nvec, 'vq_search_mag', gui_en=0);

  figure(3); clf; compare_hist("hpf150 K=64 SD Histograms", sdpf_mse, sdpf_gain, sdpf_mag)
endfunction


function short_detailed_test(train_func, test_func)
  K = 10;
  load surf_train_120;
  load surf_all;
  NtrainVec = 1000;
  NtestVec = 100;
  trainvec = surf_train_120(1:NtrainVec,1:K);
  testvec = surf_all(1:NtestVec,1:K);

  Nvec = 9;  % we can plot all vectors on one screen of subplots
  
  [idx vq]  = kmeans2(trainvec, Nvec,
                      "start", "sample", 
                      "emptyaction", "singleton", 
                      "search_func", train_func);

  sdpf = run_test(vq, testvec, Nvec, test_func, gui_en=1);
endfunction


% Some contrived examples to test VQ training

function test_training_mse
  K = 3; NtrainVec = 10; Nvec = 2;  

  trainvec = ones(NtrainVec,K);
  trainvec(NtrainVec/2+1:NtrainVec,:) = -1;
  
  [idx vq]  = kmeans2(trainvec, Nvec, "emptyaction", "singleton");

  ok = find(vq == [1 1 1]) && (find(vq == [-1 -1 -1]));
  printf("ok: %d\n", ok);
endfunction


function test_training_gain
  K = 3; NtrainVec = 10; Nvec = 2;  

  % Vectors that are the same, but offset from each other via a linear
  % term. Training algorithm should map these all to the "same"
  % vector.

  trainvec = ones(NtrainVec,K);
  for v=2:NtrainVec/2
   trainvec(v,:) += v*ones(1,1:K);
  end

  % Second set of "identical" vectors  except for gain offset

  for v=NtrainVec/2+1:NtrainVec
   trainvec(v,:) = [1 0 -1] + (v-NtrainVec/2-1)*ones(1,1:K);
  end

  [idx vq]  = kmeans2(trainvec, Nvec, 
                      "start", "sample", 
                      "emptyaction", "singleton", 
                      "search_func", "vq_search_gain");

  % check we get a vq table of two vectors that are linear offset [1 1 1] and [1 0 -1]

  tol = 0.001; ok = 0;
  for i=1:Nvec
    diff = vq(i,:) - [1 1 1];
    if std(diff) < tol, ok++; end;
    diff = vq(i,:) - [1 0 -1];
    if std(diff) < tol, ok++; end;
  end
  if ok == 2, printf("gain: OK\n"); end;
endfunction


function test_training_mag
  K = 3; NtrainVec = 10; Nvec = 2;  

  % Given a vector x, create a set of training data y = m*x + c, with x
  % modified by a magnitude and linear term.  Each vector has a
  % different mag and linear term.

  trainvec = zeros(NtrainVec,K);
  for v=1:2:NtrainVec
    trainvec(v,:) = v*[1 2 3] + 2*v;
  end

  % another set of "identical" vectors, mapped by different magnitude and linear terms,
  % alternated with the frist set so we can use the "start:first" to populate the VQ

  for v=2:2:NtrainVec
    trainvec(v,:) = cos(v)*[2 -1 2] -2*v;
  end

  trainvec

  [idx vq]  = kmeans2(trainvec, Nvec, 
                      "start", "first", 
                      "search_func", "vq_search_mag");

  vq

  % todo: how to auto test?  Need to solve same euqations?
#}
  tol = 0.001; ok = 0;
  for i=1:Nvec
    diff = vq(i,:) - [1 1 1];
    if std(diff) < tol, ok++; end;
    diff = vq(i,:) - [1 0 -1];
    if std(diff) < tol, ok++; end;
  end
  if ok == 2, printf("gain: OK\n"); end;
#}
endfunction


% ---------------------------------------------------------

format; more off;
rand('seed',1);    % kmeans using rand for initial population,
                   % we want same results on every run

%long_tests(quick_check=1)
short_detailed_test('vq_search_mag', 'vq_search_mag');
%test_training_mag



