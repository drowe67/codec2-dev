% quant_est.m
% Estimate number of bits to quantise a data set, using ideas in [1]
% [1] Makhoul, "Vector Quantisation in Speech Coding", 1985
%
% This idea is to use these estimates to:
% a) compare to VQ results of same dataset
% b) get a quick indication of how preprocessing data set will affect #bits
%    required for a target variance.  It can be argued that delta bits 
%    compared to VQ will be consistent for a given target_var.
%
% usage:
%   octave:> quant_est('../build_linux/train_120_b.f32',13,1)
 
function quant_est(fn,target_var,kmeans_en=0,K=30)
  x = load_f32(fn,K);

  % remove low energy vectors and mean
  Kst=0; Ken=29; x = reshape(x(:,Kst+1:Ken+1),length(x),Ken-Kst+1);
  K = Ken-Kst+1;
  x = x(find(mean(x,2)>10),:);
  x -= mean(x,2);

  printf("                           Bits  Var\n");
  % quantise using eigen rotation and scalar quantisers 
  lambda = eig(cov(x))';
  [mean_var nbits] = allocate_bits(lambda,target_var);
  printf("eigen rotation and scalar: %2d    %5.2f\n", sum(nbits), mean_var);
 
  % quantise using DCT and scalar quantisers 
  pkg load signal
  dct_var = var(dct2(x));
  [mean_var nbits] = allocate_bits(dct_var,target_var);
  printf("DCT and scalar...........: %2d    %5.2f\n", sum(nbits), mean_var);

  if kmeans_en
    M=512;
    [idx centers]=kmeans(x,M); 
    mean_var = mean(var(centers(idx,:) - x));
    printf("kmeans VQ................: %2d    %5.2f\n", log2(M), mean_var);
  end  
endfunction 

% Keep allocating bits to simulated scalar qauntisers until average distortion 
% reaches we reach target_var. We assume each bit drops variance by 6dB.  A
% more accurate approach would be to actually train quantisers
function [mean_var nbits] = allocate_bits(var, target_var)
  K = length(var);
  nbits = zeros(1,K); nbits_total = 0; 
  mean_var = mean(var);
  delta_var = zeros(1,K);
  while mean_var > target_var
    for k=1:K
      delta_var(k) = var(k)/4;
    end
    
    [S I] = sort(delta_var,'descend');
    nbits(I(1))++; nbits_total++;
    var(I(1)) = delta_var(I(1));
    mean_var = mean(var);
    %printf("nbits: %d mean var: %f\n", nbits_total, mean_var);
  endwhile
endfunction
