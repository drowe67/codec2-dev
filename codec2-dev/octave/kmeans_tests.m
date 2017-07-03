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

function [idx contrib errors test_ g mg] = vq_search_mse(vq, data)
  [nVec nCols] = size(vq);
  nRows = length(data);
  
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

  g = mg = 1; % dummys for this function
endfunction


% abs() search with a linear gain term

function [idx contrib errors test_ g mg] = vq_search_gain(vq, data)
  
  [nVec nCols] = size(vq);
  nRows = length(data);

  error = zeros(1,nVec);
  g = zeros(nRows, nVec);
  diff  = zeros(nVec, nCols);
  errors = zeros(1, nRows);
  idx = zeros(1, nRows);
  contrib = zeros(nRows, nCols);
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

    idx(f) = min_ind;

    errors(f) = mn; 
    contrib(f,:) = test_(f,:) = vq(min_ind,:) + g(f,min_ind);
  end
  mg = 1;
endfunction


% abs() search with a linear plus ampl scaling term

function [idx contrib errors test_ g mg] = vq_search_mag(vq, data)
  [nVec nCols] = size(vq);
  nRows = length(data);
  
  g = mg = zeros(nRows, nVec);
  diff  = zeros(nVec, nCols);
  errors = zeros(1, nRows);
  idx = error = zeros(1, nVec);
  contrib = zeros(nRows, nCols);
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

    contrib(f,:) = test_(f,:) = mg(f,min_ind) * vq(min_ind,:) + g(f,min_ind);
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

  [idx contrib errors test_ g mg] = feval(search_func, vq, test);
  
  % sd over time
 
  sd_per_frame = zeros(nRows,1);
  for i=1:nRows
    sd_per_frame(i) = mean(abs(test(i,:) - test_(i,:)));
  end

  printf("%18s mean SD: %3.2f dB\n", search_func, mean(sd_per_frame));
  
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


function sd = three_tests(sim_in)
  trainvec = sim_in.trainvec;
  testvec = sim_in.testvec;
  Nvec = sim_in.Nvec;
  train_func = sim_in.train_func;

  printf("  Nvec: %d\n", Nvec);

  [idx vq] = kmeans(trainvec, Nvec, 
                    "start", "sample", 
                    "emptyaction", "singleton",
                    "search_func", train_func);

  sd_mse  = run_test(vq, testvec, Nvec, 'vq_search_mse', gui_en=0);
  sd_gain = run_test(vq, testvec, Nvec, 'vq_search_gain', gui_en=0);
  sd_mag = run_test(vq, testvec, Nvec, 'vq_search_mag', gui_en=0);

  sd = [sd_mse sd_gain sd_mag];
endfunction

function plot_sd_results(title_str, fg, offset, sd)
  figure(fg); clf;  
  samples = offset+[1 4 7];
  bits = log2([64 128 256]);
  errorbar(bits-0.1, mean(sd(:,samples)), std(sd(:, samples),[]),'b+-; MSE search;');
  hold on;
  samples = offset+[2 5 8];
  errorbar(bits+0.0, mean(sd(:,samples)), std(sd(:,samples),[]),'g+-;Gain search;');
  samples = offset+[3 6 9];
  errorbar(bits+0.1, mean(sd(:,samples)), std(sd(:,samples),[]),'r+-;Mag search;');
  hold off;
  xlabel('VQ size (bits)')
  ylabel('mean SD (dB)');
  title(title_str);
endfunction


function plot_sd_results2(title_str, fg, sd)
  figure(fg); clf;  
  samples = 0+[3 6 9];
  bits = log2([64 128 256]);
  errorbar(bits-0.1, mean(sd(:,samples)), std(sd(:, samples),[]),'b+-; MSE train;');
  hold on;
  samples = 9+[3 6 9];
  errorbar(bits+0.0, mean(sd(:,samples)), std(sd(:,samples),[]),'g+-;Gain train;');
  samples = 18+[3 6 9];
  errorbar(bits+0.1, mean(sd(:,samples)), std(sd(:,samples),[]),'r+-;Mag train;');
  hold off;
  xlabel('VQ size (bits)')
  ylabel('mean SD (dB)');
  title(title_str);
endfunction


function sd = long_tests(quick_check=0)
  num_cores = 4;
  K = 10;
  load surf_train_120; load surf_all;
  load surf_train_120_hpf150; load surf_all_hpf150;

  if quick_check
    NtrainVec = 1000;
    NtestVec = 100;
  else
    NtrainVec = length(surf_train_120);
    NtestVec = length(surf_all);
  end

  trainvec = surf_train_120(1:NtrainVec,1:K);
  testvec = surf_all(1:NtestVec,1:K);
  trainvec_hpf150 = surf_train_120_hpf150(1:NtrainVec,1:K);
  testvec_hpf150 = surf_all_hpf150(1:NtestVec,1:K);

  sim_in.trainvec = trainvec; sim_in.testvec = testvec; sim_in.Nvec = 64; sim_in.train_func = 'vq_search_mse';
  sim_in_vec(1:3) = sim_in;
  sim_in_vec(2).Nvec = 128; sim_in_vec(3).Nvec = 256;
  sd =  pararrayfun(num_cores, @three_tests, sim_in_vec);
  plot_sd_results("MSE Training", 1, 0, sd)

  for i=1:3, sim_in_vec(i).trainvec = trainvec_hpf150; sim_in_vec(i).testvec = testvec_hpf150; end;
  sd = [sd pararrayfun(num_cores, @three_tests, sim_in_vec)];
  plot_sd_results("MSE training 150Hz HPF", 2, 9, sd)

  for i=1:3, sim_in_vec(i).train_func = 'vq_search_gain'; end;
  sd = [sd pararrayfun(num_cores, @three_tests, sim_in_vec)];
  plot_sd_results("Gain training 150Hz HPF", 3, 9, sd)

  for i=1:3, sim_in_vec(i).train_func = 'vq_search_mse'; end;
  sd = [sd pararrayfun(num_cores, @three_tests, sim_in_vec)];
  plot_sd_results("Mag training 150Hz HPF", 4, 9, sd)

  plot_sd_results2("Mag search 150Hz HPF", 5, sd)
  
  figure(6); clf; compare_hist("Mag train Nvec=64", sd(:,3), sd(:,12), sd(:,21));
endfunction


function short_detailed_test(train_func, test_func)
  K = 10;
  load surf_train_120_hpf150;
  load surf_all;
  NtrainVec = 1000;
  NtestVec = 100;
  trainvec = surf_train_120_hpf150(1:NtrainVec,1:K);
  testvec = surf_all(1:NtestVec,1:K);

  Nvec = 64;  % we can plot all vectors on one screen of subplots
  
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
  trainvec(2:2:NtrainVec,:) = -1;
  
  [idx vq]  = kmeans2(trainvec, Nvec,
                      "start", "first", 
                      "emptyaction", "singleton",
                      "search_func", "vq_search_mse");

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

%short_detailed_test('vq_search_mag', 'vq_search_mag');
sd = long_tests(quick_check=0);
%test_training_mag



