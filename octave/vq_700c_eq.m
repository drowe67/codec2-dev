% vq_700c.m
% David Rowe May 2019
%
% Researching Codec 2 700C VQ equaliser ideas

melvq;

% general purpose plot function for looking at averages of K-band
% sequences in scripts dir and VQs:
%   vq_700c_plots({"hts2a.f32" "vk5qi.f32" "train_120_1.txt"})

function vq_700c_plots(fn_array)
  nb_features = 41
  K = 20
  figure(1); clf; hold on; axis([1 20 -20 40]); title('Max Hold');
  figure(2); clf; hold on; axis([1 20 -20 30]); title('Average'); 
  for i=1:length(fn_array)
    [dir name ext] = fileparts(fn_array{i});
    if strcmp(ext, ".f32")
      % f32 feature file
      fn = sprintf("../script/%s_feat%s", name, ext)
      feat = load_f32(fn , nb_features);
      bands = feat(:,2:K+1);
    else
      % text file (e.g. existing VQ)
      bands = load(fn_array{i});
    end
    figure(1); plot(max(bands),'linewidth', 5);
    figure(2); plot(mean(bands),'linewidth', 5);
  end
  figure(1); legend(fn_array);
  figure(2); legend(fn_array);
endfunction


% single stage vq a target matrix

function errors = vq_targets(vq, targets)
  errors = [];
  for i=1:length(targets)
    [mse_list index_list] = search_vq(vq, targets(i,:), 1);
    error = targets(i,:) - vq(index_list(1),:);
    errors = [errors; error];
  end
endfunction


% two stage mbest vq a target matrix

function [errors targets_] = vq_targets2(vq1, vq2, targets)
  vqset(:,:,1)= vq1; vqset(:,:,2)=vq2; m=5;
  [errors targets_] = mbest(vqset, targets, m);
endfunction


% given target and vq matrices, estimate eq via two metrics

function [eq1 eq2] = est_eq(vq, targets)
  [ntargets K] = size(targets);
  [nvq K] = size(vq);
  
  eq1 = zeros(1,K);  eq2 = zeros(1,K);
  for i=1:length(targets)
    [mse_list index_list] = search_vq(vq, targets(i,:), 1);

    % eq metric 1: average of error for best VQ entry
    eq1 += targets(i,:) - vq(index_list(1),:);
    
    % eq metric 2: average of error across all VQ entries
    for j=1:nvq
      eq2 += targets(i,:) - vq(j,:);
    end
  end

  eq1 /= ntargets;
  eq2 /= (ntargets*nvq);
endfunction

function save_f32(fn, m)
  f=fopen(fn,"wb");
  [r c] = size(m);
  mlinear = reshape(m', 1, r*c);
  fwrite(f, mlinear, 'float32');
  fclose(f);
endfunction

function [targets e] = load_targets(fn_target_f32)
  nb_features = 41;
  K = 20;

  % .f32 files are in scripts directory, first K values rate_K_no_mean vectors
  [dir name ext] = fileparts(fn_target_f32);
  fn = sprintf("../script/%s_feat.f32", name);
  feat = load_f32(fn, nb_features);
  e = feat(:,1);
  targets = feat(:,2:K+1);
endfunction


function table_across_samples
  K = 20;

  % VQ is in .txt file in this directory, we have two to choose from.  train_120 is the Codec 2 700C VQ,
  % train_all_speech was trained up from a different, longer database, as a later exercise
  #vq_name = "train_120";
  vq_name = "train_all_speech";  
  vq1 = load(sprintf("%s_1.txt", vq_name));
  vq2 = load(sprintf("%s_2.txt", vq_name));
  
  printf("-------------------------------------------------------------------\n");
  printf("Sample            Initial  vqstg1  vqstg1_eq   vqsgt2  vq2stg_eq\n");
  printf("-------------------------------------------------------------------\n");
            
  fn_targets = {"hts1a" "hts2a" "cq_ref" "ve9qrp_10s" "vk5qi" "c01_01_8k" "ma01_01" "cq_freedv_8k"};
  %fn_targets = {"hts1a"};
  for i=1:length(fn_targets)

    % load target and estimate eq
    [targets e] = load_targets(fn_targets{i});
    eq = est_eq(vq1, targets);

    % first stage VQ
    errors1 = vq_targets(vq1, targets);
    errors1_eq = vq_targets(vq1, targets-eq);
    % two stage mbest VQ
    [errors2 targets_] = vq_targets2(vq1, vq2, targets);
    [errors2_eq targets_eq_] = vq_targets2(vq1, vq2, targets-eq);

    % save to .f32 files for listening tests
    if strcmp(vq_name,"train_120")
      save_f32(sprintf("../script/%s_vq2.f32", fn_targets{i}), targets_);
      save_f32(sprintf("../script/%s_vq2_eq.f32", fn_targets{i}), targets_eq_);
    else
      save_f32(sprintf("../script/%s_vq2_as.f32", fn_targets{i}), targets_);
      save_f32(sprintf("../script/%s_vq2_as_eq.f32", fn_targets{i}), targets_eq_);
    end 
    printf("%-17s %6.2f  %6.2f  %6.2f     %6.2f  %6.2f\n", fn_targets{i},
            var(targets(:)), var(errors1(:)), var(errors1_eq(:)),
            var(errors2(:)), var(errors2_eq(:)));
   end
endfunction


% interactve, menu driven frame by frame plots

function interactive(fn_vq_txt, fn_target_f32)
  K = 20;
  vq = load("train_120_1.txt");
  [targets e] = load_targets(fn_target_f32);
  hi_energy_frames = find(e>20);
  [eq1 eq2] = est_eq(vq, targets);
  [eq1_hi eq2_hi] = est_eq(vq, targets(hi_energy_frames,:));
  
  figure(1); clf;
  mesh(e+targets)
  figure(2); clf;
  plot(eq1,'b;eq1;')
  hold on;
  plot(eq2,'g;eq2;');
  plot(eq1_hi,'r;eq1 hi;');
  plot(eq2_hi,'c;eq2 hi;');
  hold off;
  figure(3); clf; plot(e(:,1))

  % enter single step loop
  f = 20; neq = 0; eq=zeros(1,K);
  do 
    figure(4); clf;
    t = targets(f,:) - eq;
    [mse_list index_list] = search_vq(vq, t, 1);
    error = t - vq(index_list(1),:);
    plot(e(f)+t,'b;target;');
    hold on;
    plot(e(f)+vq(index_list,:),'g;vq;');
    plot(error,'r;error;');
    plot([1 K],[e(f) e(f)],'--')
    hold off;
    axis([1 K -20 80])
    % interactive menu 

    printf("\r f: %2d eq: %d ind: %3d var: %3.1f menu: n-next  b-back  e-eq q-quit", f, neq, index_list(1), var(error));
    fflush(stdout);
    k = kbhit();

    if k == 'n' f+=1; end
    if k == 'e'
      neq++;
      if neq == 3 neq = 0; end
      if neq == 0 eq = zeros(1,K); end
      if neq == 1 eq = eq1; end
      if neq == 2 eq = eq2; end
    end
    if k == 'b' f-=1; end
  until (k == 'q')
  printf("\n");
endfunction

more off
%interactive("train_120_1.txt", "ve9qrp_10s.f32")
table_across_samples;
%vq_700c_plots({"hts1a.f32" "hts2a.f32" "ve9qrp_10s.f32" "ma01_01.f32" "train_120_1.txt"})
