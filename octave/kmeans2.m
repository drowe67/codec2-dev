## Copyright (C) 2011 Soren Hauberg
## Copyright (C) 2012 Daniel Ward
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see .
## -*- texinfo -*-
## @deftypefn {Function File} {[@var{idx}, @var{centers}] =} kmeans2 (@var{data}, @var{k}, @var{param1}, @var{value1}, @dots{})
## K-means clustering.
##
## @seealso{linkage}
## @end deftypefn

function [classes, centers, sumd, D] = kmeans2 (data, k, varargin)

  [reg, prop] = parseparams (varargin);

  ## defaults for options

  emptyaction = "error";
  start       = "sample";

  # used for getting the number of samples

  nRows = rows (data);

  ## used to hold the distances from each sample to each class

  D = zeros (nRows, k);

  # used for convergence of the centroids

  err = 1;

  # initial sum of distances

  sumd = Inf;

  # default search function, can  be over-ridden by user supplied function

  search_func = @vq_search_mse;

  ## Input checking, validate the matrix and k

  if (!isnumeric (data) || !ismatrix (data) || !isreal (data))
    error ("kmeans: first input argument must be a DxN real data matrix");
  elseif (!isscalar (k))
    error ("kmeans: second input argument must be a scalar");
  endif

  if (length (varargin) > 0)

    ## check for the ‘emptyaction’ property

    found = find (strcmpi (prop, "emptyaction") == 1);
    if found 
      switch (lower (prop{found+1}))
        case "singleton"
          emptyaction = "singleton";
        otherwise
          error ("kmeans: unsupported empty cluster action parameter");
      endswitch
    end

    ## check for the ‘search_func’ property, user defined vq_search function

    found = find (strcmpi (prop, "search_func") == 1);
    if found
      search_func = prop{found+1};
    end

    ## check for the ‘start’ property

    found = find (strcmpi (prop, "start") == 1);
    if found
      switch (lower (prop{found+1}))
        case "sample"
          idx = randperm (nRows) (1:k);
          centers = data (idx, :);
        case "first"
          centers = data (1:k, :);
        otherwise
          error ("kmeans: unsupported initial clustering parameter");
      endswitch
    end
  end

  ## Run the algorithm

  while err > .001
    classes = feval(search_func, centers, data);

    ## Calculate new centroids

    for i = 1:k

      ## Get binary vector indicating membership in cluster i

      membership = (classes == i);

      ## Check for empty clusters

      if (sum (membership) == 0)
        switch emptyaction

          ## if ‘singleton’, then find the point that is the
          ## farthest and add it to the empty cluster

          case 'singleton'
           idx=maxCostSampleIndex (data, centers(i,:));
           classes(idx) = i;
           membership(idx)=1;
         ## if ‘error’ then throw the error
          otherwise
            error ("kmeans: empty cluster created");
        endswitch
      endif ## end check for empty clusters

      ## update the centroids

      members = data(membership, :);
      centers(i, :) = sum(members,1)/size(members,1);
    endfor

    ## calculate the difference in the sum of distances

    err  = sumd - objCost (data, classes, centers);

    ## update the current sum of distances

    sumd = objCost (data, classes, centers);

  endwhile
endfunction


function idx = vq_search_mse(vq, data)
    [nVec nCols] = size(vq);
    nRows = length(data);

    error = zeros(1,nVec);
    idx = zeros(1, nRows);

    for f=1:nRows
      target = data(f,:);
      for i=1:nVec
        diff = target - vq(i,:);
        error(i) = diff * diff';
      end
      [mn min_ind] = min(error);
      idx(f) = min_ind;
     end
endfunction


## calculate the sum of distances

function obj = objCost (data, classes, centers)
  obj = 0;
    for i=1:rows (data)
      obj = obj + sumsq (data(i,:) - centers(classes(i),:));
    endfor
endfunction

function idx = maxCostSampleIndex (data, centers)
  cost = 0;
  for idx = 1:rows (data)
    if cost < sumsq (data(idx,:) - centers)
      cost = sumsq (data(idx,:) - centers);
    endif
  endfor
endfunction

%!demo
%! ## Generate a two-cluster problem
%! C1 = randn (100, 2) + 1;
%! C2 = randn (100, 2) – 1;
%! data = [C1; C2];
%!
%! ## Perform clustering
%! [idx, centers] = kmeans (data, 2);
%!
%! ## Plot the result
%! figure
%! plot (data (idx==1, 1), data (idx==1, 2), ‘ro’);
%! hold on
%! plot (data (idx==2, 1), data (idx==2, 2), ‘bs’);
%! plot (centers (:, 1), centers (:, 2), ‘kv’, ‘markersize’, 10);
%! hold off
