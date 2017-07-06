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

%----------------------------------------------------------------------
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
  mg = 1; sl = 0;
endfunction


%----------------------------------------------------------------------
% abs() search with a linear plus ampl scaling term

function [idx contrib errors test_ g mg sl] = vq_search_mag(vq, data)
  [nVec nCols] = size(vq);
  nRows = rows(data);
  
  g = mg = zeros(nRows, nVec);
  diff  = zeros(nVec, nCols);
  errors = zeros(1, nRows);
  idx = error = zeros(1, nVec);
  contrib = zeros(nRows, nCols);
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
  end
  sl = 0;
endfunction


%----------------------------------------------------------------------
% abs() search with a linear, ampl scaling, and slope term

function [idx contrib errors test_ g mg sl] = vq_search_slope(vq, data)
  [nVec nCols] = size(vq);
  nRows = rows(data);
  
  g = mg = sl = zeros(nRows, nVec);
  diff  = zeros(nVec, nCols);
  errors = zeros(1, nRows);
  idx = error = zeros(1, nVec);
  contrib = zeros(nRows, nCols);
  test_ = zeros(nVec, nCols);

  weights = ones(1,nCols);

  for f=1:nRows
    target = data(f,:);
    %target = 2*vq(1,:)+1;

    for i=1:nVec
      % work out gain, amp and slope for best match, 3 unknowns, 3 equations

      A = [sum(vq(i,:))      nCols        sum(1:nCols); ...
          vq(i,:)*vq(i,:)'   sum(vq(i,:)) (1:nCols)*vq(i,:)'; ...
          (1:nCols)*vq(i,:)' sum(1:nCols) (1:nCols)*(1:nCols)'];
      c = [sum(target) target*vq(i,:)' target*(1:nCols)']';
      b = inv(A)*c;

      mg(f,i) = b(1); g(f,i) = b(2); sl(f,i) = b(3);
      diff(i,:) = target - (mg(f,i)*vq(i,:) + g(f,i) + sl(f,i)*(1:nCols));
      diff(i,:) .* weights;

      % abs in dB is MSE in linear

      error(i) = mean(abs(diff(i,:)));

      %printf("f: %d i: %d mg: %f g: %f sl: %f error: %f\n", f, i, mg(f,i), g(f,i), sl(f,i), error(i));
    end
    
    [mn min_ind] = min(error);
    errors(f) = mn; 
    idx(f) = min_ind;

    contrib(f,:) = test_(f,:) = mg(f,min_ind)*vq(min_ind,:) + g(f,min_ind) + sl(f,min_ind)*(1:nCols);
  end

endfunction


%----------------------------------------------------------------------
%
% Functions to support simulation of different VQ training and testing
%
%----------------------------------------------------------------------

% evaluate database test using vq, with selectable search function.  Can be operated in
% GUI mode to analyse in fine detail or batch mode to evaluate lots of data.

function sd_per_frame = run_test(vq, test, nVec, search_func, gui_en = 0)

  % Test VQ using test data -----------------------

  % Each test must output:
  %   errors: - error for best vector
  %   idx...: - index of best vector
  %   test_.: - best approximation of target test database
  
  [nRows nCols] = size(test);

  [idx contrib errors test_ g mg sl] = feval(search_func, vq, test);

  % sd over time
 
  sd_per_frame = zeros(nRows,1);
  for i=1:nRows
    sd_per_frame(i) = mean(abs(test(i,:) - test_(i,:)));
  end

  printf("average sd: %3.2f\n", mean(sd_per_frame));
  %printf("%18s nVec: %d SD: %3.2f dB\n", search_func, nVec, mean(sd_per_frame));
  
  % Optional GUI with U to anlayse VQ perf ---------------------------------------------------

  if gui_en
    gui_mode = "sort";     % "sort" or "fbf"
    sort_order = "worst";  % "worst" or "best"
    f = 1;

    % plots sd over time and histogram

    figure(1); clf; subplot(211); plot(sd_per_frame); title('SD'); 
    subplot(212); [yy xx] = hist(sd_per_frame); plot(xx,yy,'+-');
  end

  while gui_en

    % display m frames, printing some stats, plotting vector to give visual idea of match

    figure(2); clf;

    % try this in 'ascend' or 'descend' mode to best or worst frames

    if strcmp(gui_mode, "sort")
      if strcmp(sort_order, "best")
        [errors_dec frame_dec] = sort(errors, "ascend");
      end
      if strcmp(sort_order, "worst")
        [errors_dec frame_dec] = sort(errors, "descend");
      end
       m = 4;
    else
      m = 1;
      errors_dec = errors(f); frame_dec(1) = f;
    end

    for i=1:m

      % build up VQ legend

      af = frame_dec(i); aind = idx(af);

      ag = 0; amg = 1; asl = 0;
      if strcmp(search_func, "vq_search_gain")
        ag = g(af,aind);
        l3 = sprintf("g: %3.2f", ag);
      end
      if strcmp(search_func, "vq_search_mag")
        ag = g(af,aind);
        amg = mg(af,aind);
        l3  = sprintf("g: %3.2f mg: %3.2f", ag, amg);
      end
      if strcmp(search_func, "vq_search_slope")
        ag = g(af,aind);
        amg = mg(af,aind);
        asl = sl(af,aind);
        l3 = sprintf("g: %3.2f mg: %3.2f sl: %3.2f", ag, amg, asl);
      end

      % plot target 

      subplot(sqrt(m), sqrt(m),i);
      l1 = sprintf("b-;fr %d sd: %3.2f;", af, errors_dec(m));
      plot(test(af,:), l1); 

      % plot vq vector and modified version

      hold on; 
      l2 = sprintf("g-+;ind %d;", aind); 
      plot(vq(aind, :), l2); 
      l3 = sprintf("g-o;%s;",l3);
      plot(amg*vq(aind, :) + ag + asl*(1:nCols), l3); 
      hold off;
      axis([1 nCols -10 40]);
    end

    % interactive menu ------------------------------------------

    if strcmp(gui_mode, "sort")
      printf("\rmenu: m-mode[%s] o-order[%s] q-quit", gui_mode, sort_order);
    end
    if strcmp(gui_mode, "fbf")
      printf("\rmenu: m-mode[%s] frame: %d  n-next  b-back  q-quit", gui_mode, f);
    end
    fflush(stdout);
    k = kbhit();
    
    if k == 'm'
      if strcmp(gui_mode, "sort")
        gui_mode = "fbf"; f = frame_dec(1); 
      else
        gui_mode = "sort";
      end;
    end
    if k == 'o'
       if strcmp(sort_order, "worst")
         sort_order = "best"; 
       else
         sort_order = "worst";
      end;
    end
    if k == 'n', f = f + 1; endif;
    if k == 'b', f = f - 1; endif;
    if k == 'q', gui_en = 0; printf("\n"); endif;
  end  
  
endfunction


function search_func = get_search_func(short_name)
  if strcmp(short_name, "mse")
    search_func = 'vq_search_mse';
  end
  if strcmp(short_name, "gain")
    search_func = 'vq_search_gain';
  end
  if strcmp(short_name, "mag")
    search_func = 'vq_search_mag';
  end
  if strcmp(short_name, "slope")
    search_func = 'vq_search_slope';
  end
end


%----------------------------------------------------------------------
% Train up a  VQ and run one or mode tests

function [sd des] = train_vq_and_run_tests(sim_in)
  Nvec = sim_in.Nvec;

  train_func = get_search_func(sim_in.train_func_short);
  printf("    Nvec %d kmeans start\n", Nvec);
  [idx vq] = feval(sim_in.kmeans,
                   sim_in.trainvec, Nvec, 
                   "start", "sample", 
                   "emptyaction", "singleton",
                   "search_func", train_func);
  printf("    Nvec %d kmeans end\n", Nvec);

  sd = []; des = []; tests = sim_in.tests;
  for t=1:length(cellstr(tests))
    test_func_short = char(cellstr(tests)(t));
    test_func = get_search_func(test_func_short);
    asd = run_test(vq, sim_in.testvec, Nvec, test_func); 
    sd = [sd asd];
    ades.Nvec = Nvec; 
    ades.train_func_short = sim_in.train_func_short;
    ades.test_func_short  = test_func_short;
    ades_str = sprintf("Nvec: %3d train: %-5s test: %-5s", 
                       ades.Nvec, ades.train_func_short, ades.test_func_short); 
    des = [des ades];
    printf("    %s SD: %3.2f\n", ades_str, mean(asd));
  end

endfunction

% Search for test that matchs train/test and any Nvec and constructs points for error bar plot

function [y x leg] = search_tests(sd, desc, train, test)
  nTests = length(desc)
  x = []; y = []; leg = [];
  for i=1:nTests
    if strcmp(desc(i).train_func_short, train) && strcmp(desc(i).test_func_short, test)
      printf("i: %2d Nvec: %3d train: %5s test: %5s\n", i, desc(i).Nvec, desc(i).train_func_short, desc(i).test_func_short);
      x = [x; log2(desc(i).Nvec)];
      y = [y; mean(sd(:,i)) std(sd(:,i))];
    end
  end
endfunction


function plot_sd_results(fg, title_str, sd, desc, train_list, test_list)
  figure(fg); clf;  

  nlines = 0; inc = -0.1;
  for i=1:length(cellstr(train_list))
    train_func_short = char(cellstr(train_list)(i));
    for j=1:length(cellstr(test_list))
      test_func_short = char(cellstr(test_list)(j));
      [y x] = search_tests(sd, desc, train_func_short, test_func_short);
      leg = sprintf("o-%d;train: %5s test: %5s;", nlines, train_func_short, test_func_short);
      if nlines, hold on; end;
      x += inc; inc += 0.1;            % separate x coords a bit to make errors bars legible
      errorbar(x, y(:,1), y(:,2), leg);
      nlines++;
    end
    if nlines>1, hold off; end;
  end

  xlabel('VQ size (bits)')
  ylabel('mean SD (dB)');
  title(title_str);
endfunction


% Search for test that matches all fields for histogram plots

function testNum = search_tests_Nvec(sd, desc, train, test, Nvec)
  nTests = length(desc);
  testNum = 0;
  for i=1:nTests
    if strcmp(desc(i).train_func_short, train) && strcmp(desc(i).test_func_short, test) && desc(i).Nvec == Nvec
      printf("i: %2d Nvec: %3d train: %5s test: %5s\n", i, desc(i).Nvec, desc(i).train_func_short, desc(i).test_func_short);
      testNum = i;
    end
  end
endfunction

%----------------------------------------------------------------------
%
% Plot histograms of SDs for comparison. Each col of sd_per_frame
% has the results of one test, number of cols % is number of tests.  leg
% is a col vector with one legend string for each test.

function compare_hist(fg, atitle, sd, desc)
  figure(fg); clf;
  [nRows nCols] = size(sd);
  for c=1:nCols
    [yy, xx] = hist(sd(:, c));
    if c == 2, hold on; end;
    leg = sprintf("o-%d;train: %5s test: %5s;", c-1, desc(c).train_func_short, desc(c).test_func_short);
    plot(xx, yy, leg);
  end
  if nCols > 1, hold off; end;
  title(atitle)
end

%----------------------------------------------------------------------
% Run a bunch of long tests in parallel and plot results
%

% This test tries a bunch of training and test search combinations on
% the first 1kHz to see if there is any advantage in using the higher
% order search techniques in the VQ training.  So far it appears
% not....  However the higher order search routines do give great
% results compared to a basic MSE (or oder) search.

function [sd desc] = long_tests_1_10(quick_check=0)
  num_cores = 4;
  K = 10;
  load surf_train_120_hpf150; load surf_all_hpf150;

  if quick_check
    NtrainVec = 1000;
    NtestVec = 100;
  else
    NtrainVec = length(surf_train_120_hpf150);
    NtestVec = length(surf_all_hpf150);
  end

  trainvec = surf_train_120_hpf150(1:NtrainVec,1:K);
  testvec = surf_all_hpf150(1:NtestVec,1:K);

  % build up a big array of tests to run ------------------------

  sim_in.trainvec = trainvec; sim_in.testvec = testvec;
  sim_in.kmeans = "kmeans2";  % slow but supports different search functions
  sim_in.Nvec = 64; 
  sim_in.train_func_short = "mse"; 
  sim_in.tests = ["mse"; "gain"; "mag"; "slope"];

  % Test1: mse training, 64, 128, 256

  sim_in_vec(1:3) = sim_in;
  sim_in_vec(2).Nvec = 128; sim_in_vec(3).Nvec = 256;

  test_list = sim_in_vec;

  % Test2: gain training, 64, 128, 256

  for i=1:3, sim_in_vec(i).train_func_short = "gain"; end;
  test_list = [test_list sim_in_vec];

  % Test3: mag training, 64, 128, 256

  for i=1:3, sim_in_vec(i).train_func_short = "mag"; end;
  test_list = [test_list sim_in_vec];

  % Test4: slope training, 64, 128, 256

  for i=1:3, sim_in_vec(i).train_func_short = "slope"; end;
  test_list = [test_list sim_in_vec];

  % run test list in parallel

  [sd desc] =  pararrayfun(num_cores, @train_vq_and_run_tests, test_list);
  %train_vq_and_run_tests(test_list(1));
  
  % Plot results -----------------------------------------------

  fg = 1;
  plot_sd_results(fg++, "MSE Training", sd, desc, "mse", ["mse"; "gain"; "mag"; "slope"]);
  plot_sd_results(fg++, "Gain Training", sd, desc, "gain", ["mse"; "gain"; "mag"; "slope"]);
  plot_sd_results(fg++, "Mag Training", sd, desc, "mag", ["mse"; "gain"; "mag"; "slope"]);
  plot_sd_results(fg++, "Slope Training", sd, desc, "slope", ["mse"; "gain"; "mag"; "slope"]);
  plot_sd_results(fg++, "Slope Searching", sd, desc, ["mse"; "gain"; "mag"; "slope"], "slope");

  % histogram of results from Nvec=64, all training methods, slope search

  testnum = search_tests_Nvec(sd, desc, "mse", "slope", 64);
  hist_sd = sd(:,testnum); hist_desc = desc(testnum);
  testnum = search_tests_Nvec(sd, desc, "gain", "slope", 64);
  hist_sd = [hist_sd sd(:,testnum)]; hist_desc = [hist_desc desc(testnum)];
  testnum = search_tests_Nvec(sd, desc, "mag", "slope", 64);
  hist_sd = [hist_sd sd(:,testnum)]; hist_desc = [hist_desc desc(testnum)];
  testnum = search_tests_Nvec(sd, desc, "slope", "slope", 64);
  hist_sd = [hist_sd sd(:,testnum)]; hist_desc = [hist_desc desc(testnum)];
  compare_hist(fg, "Histogram of SDs Nvec=64", hist_sd, hist_desc);

endfunction


%----------------------------------------------------------------------
% This test tries some test combinations on the 1000 to 3000 Hz range.
%
% Based on our results with long_test_1_10() we just use MSE training

function [sd desc] = long_tests_11_40(quick_check=0)
  num_cores = 4;
  K = 20; st = 11; en = 40;
  load surf_train_120_hpf150; load surf_all_hpf150;

  if quick_check
    NtrainVec = 1000;
    NtestVec = 100;
  else
    NtrainVec = length(surf_train_120_hpf150);
    NtestVec = length(surf_all_hpf150);
  end

  trainvec = surf_train_120_hpf150(1:NtrainVec,st:en);
  testvec = surf_all_hpf150(1:NtestVec,st:en);

  % build up an array of tests to run ------------------------

  sim_in.kmeans = "kmeans";
  sim_in.trainvec = trainvec; sim_in.testvec = testvec; 
  sim_in.Nvec = 64; 
  sim_in.train_func_short = "mse"; 
  sim_in.tests = ["slope"];

  % Test1: mse training, 64, 128, 256, 512

  sim_in_vec(1:4) = sim_in;
  sim_in_vec(2).Nvec = 128; sim_in_vec(3).Nvec = 256; sim_in_vec(4).Nvec = 512;

  test_list = sim_in_vec;

  % run test list in parallel

  [sd desc] =  pararrayfun(num_cores, @train_vq_and_run_tests, test_list);
  %[sd desc] =  train_vq_and_run_tests(sim_in_vec(4));

  % Plot results -----------------------------------------------

  fg = 2;
  plot_sd_results(fg++, "MSE Training 10..30", sd, desc, "mse", "slope");

#{
  % histogram of results from Nvec=64, all training methods, slope search

  testnum = search_tests_Nvec(sd, desc, "mse", "slope", 64);
  hist_sd = sd(:,testnum); hist_desc = desc(testnum);
  testnum = search_tests_Nvec(sd, desc, "gain", "slope", 64);
  hist_sd = [hist_sd sd(:,testnum)]; hist_desc = [hist_desc desc(testnum)];
  testnum = search_tests_Nvec(sd, desc, "mag", "slope", 64);
  hist_sd = [hist_sd sd(:,testnum)]; hist_desc = [hist_desc desc(testnum)];
  testnum = search_tests_Nvec(sd, desc, "slope", "slope", 64);
  hist_sd = [hist_sd sd(:,testnum)]; hist_desc = [hist_desc desc(testnum)];
  compare_hist(fg, "Histogram of SDs Nvec=64", hist_sd, hist_desc);
#}
endfunction


function vq = detailed_test(Nvec=64, st=1, en=10, quick_check=1, test_func = "vq_search_slope")
  K = en-st+1;
  load surf_train_120_hpf150;
  load surf_all_hpf150;
  if quick_check
    NtrainVec = 1000;
    NtestVec = 100;
  else
    NtrainVec = length(surf_train_120_hpf150);
    NtestVec = length(surf_all_hpf150);
    NtestVec = 100;
  end

  trainvec = surf_train_120_hpf150(1:NtrainVec,st:en);
  testvec = surf_all_hpf150(1:NtestVec,st:en);

  [idx vq]  = kmeans(trainvec, Nvec,
                     "start", "sample", 
                     "emptyaction", "singleton");

  run_test(vq, testvec, Nvec, test_func, gui_en=1);
endfunction


% ------------------------------------------------------------------
%
% Following test_* functions have some contrived examples to test VQ
% training with each training method.  However the longtest_1_10
% results suggest these training methds don't provide any improvement
% over standard MSE kmeans.

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


function test_training_slope
  K = 40; NtrainVec = 5; Nvec = 1;  

  % Given a vector x, create a set of training data y = m*x + c + sl, where:
  %  x is a vector
  %  c is a constant offset
  %  m is a magnitude multiplier on all elements of x
  %  sl is a linear slope vector

  load surf_all;
  prototype = surf_all(73,:);
  trainvec = zeros(NtrainVec,K);
  for v=1:NtrainVec
    trainvec(v,:) = 2*prototype + 1 + v*(1:K);
   end
  figure(1); clf; plot(trainvec');

  [idx vq]  = kmeans2(trainvec, Nvec, 
                      "start", "first", 
                      "search_func", "vq_search_slope");

  
  [idx contrib errors test_ g mg sl] = vq_search_slope(prototype, trainvec)
  
  % todo: how to auto test?  Need to solve same euqations?
#{
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

% choose one of these to run

detailed_test(64, 1, 20, quick_check=1);
%[sd desc] = long_tests_1_10(quick_check=1);
%[sd desc] = long_tests_11_40(quick_check=0);
%test_training_slope



