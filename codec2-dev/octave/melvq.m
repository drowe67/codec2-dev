% melvq.m
% David Rowe Aug 2015
%
% Experimenting with VQ design for mel LSPs

% todo:
% [ ] Sorting inside search what happens if we get a order issue, fix and calc SD
% [ ] Search to avoid a worst case error rather than average?  This could be included
%     in training.
% [ ] nested/staged search

1;

% train up multi-stage VQ

function vq = trainvq(training_data, Nvec, stages)

  vq = [];
  for i=1:stages
    [idx centers] = kmeans(training_data, Nvec);
    quant_error = centers(idx,:) - training_data;
    printf("mse stage %d: %f\n", i, mean(std(quant_error)));
    training_data = quant_error;
    vq(:,:,i) = centers;
  end

end


function [mse_list index_list] = search_vq(vq, target, m)

  [Nvec order] = size(vq);

  mse = zeros(1, Nvec);
 
  % find mse for each vector

  for i=1:Nvec
     mse(i) = sum((target - vq(i,:)) .^2);
  end

  % sort and keep top m matches

  [mse_list index_list ] = sort(mse);

  mse_list = mse_list(1:m);
  index_list = index_list(1:m);

endfunction


% Search multi-stage VQ, retaining m best candidates at each stage

function [res ind] = mbest(vqset, input_vecs, m)

  [Nvec order stages] = size(vqset);
  [Ninput tmp] = size(input_vecs);

  res = [];   % residual error after VQ
  ind = [];   % index of vqs

  for i=1:Ninput
  
    % first stage, find mbest candidates

    [mse_list index_list] = search_vq(vqset(:,:,1), input_vecs(i,:), m);
    cand_list = [mse_list' index_list'];
    cand_list = sortrows(cand_list,1);

    % subsequent stages ...........

    for s=2:stages

      % compute m targets for next stage, and update path

      prev_indexes = zeros(m,s-1);
      for t=1:m
        target(t,:) = input_vecs(i,:);
        for v=1:s-1
          target(t,:) -= vqset(cand_list(t,v+1),:,v);
        end
        prev_indexes(t,:) = cand_list(t,2:s);
      end
      
      % search stage s using m targets from stage s-1
      % with m targets, we do m searches which return the m best possibilities
      % so we get a matrix with one row per candidate, m*m rows total
      % prev_indexes provides us with the path through the VQs for each candidate row

      avq = vqset(:,:,s);
      cand_list = [];
      for t=1:m
        [mse_list index_list] = search_vq(avq, target(t,:), m);
        x = ones(m,1)*prev_indexes(t,:);
        cand_row = [mse_list' x index_list'];
        cand_list = [cand_list; cand_row];
      end

      % sort into m best rows
     
      cand_list = sortrows(cand_list,1);
      cand_list = cand_list(1:m,:);

    end

    % final residual
    target(1,:) = input_vecs(i,:);
    for v=1:s
      target(1,:) -= vqset(cand_list(1,v+1),:,v);
    end
    res = [res; target(1,:)];
    ind = [ind; cand_list(1,2:1+stages)];
  end

endfunction

more off;
load "../build_linux/src/all_mel.txt"
load vq;
[res1 ind] = mbest(vq, all_mel(1,:),3);
mean(std(res1))

% save text file of "vq quantised mels"
% load back into c2sim at run time
% train on continuous mels
% sorting/stability
% see how it sounds
% Goal is to get VQ sounding OK, similar to UQ at 20 or 40ms dec,
% sig better than current 700
% analysis of data, outliers, different training
