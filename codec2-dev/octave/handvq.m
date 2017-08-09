% handvq.m
% David Rowe 21 July 2017

1;


function vq = vqlow
  prototype = [-4 12 23 24 18 12 8 5 2 0];
  vq = zeros(8,10);

  vq(1,:) = [-4 23 24 18 12  8  5  2  0  0];
  vq(2,:) = [-4 12 23 24 18 12  8  5  2  0]; 
  vq(3,:) = [-4  8 18 23 24 18 12  8  5  2]; 
  vq(4,:) = [-4  6 16 21 23 24 18 12  8  5]; 
  vq(5,:) = [-4  8 23 24 18 12 24 16  4  0]; 
  vq(6,:) = [-4  8 23 24 18 12 16 24 16 12]; 
  vq(7,:) = [-4 12 12 24 12  4  4  4  4  4]; 

  figure(1); clf; plot(vq')
endfunction


% bunch of parabolic curves from traring data

function p = para(np,K=10)
  load ../build_linux/src/all_hpf150_b_log.txt
  
  % scale so MSE contribution is the same in VQ training
  
  x=std(b_log(:,2)); y=std(b_log(:,3));
  vec = [b_log(:,2)/x b_log(:,3)/y];
  [idx cent]=kmeans(vec, np);
  cent1=[cent(:,1)*x cent(:,2)*y];

  k = 1:K; k2 = k.^2;
  p = cent1*[k2; k];
endfunction


% combine vqlow with a bunch of parabolic curves

function vq = vqlow_para(np)

  % generate two sets of vectors
  
  vql = vqlow;
  [nl tmp] = size(vql);
  p   = para(np);

  % now combine all combinations

  vq = [];
  for i=1:nl
    for j=1:np
      v = vql(i,:) + p(j,:);
      vq = [vq; v];
    end
  end
endfunction


function vq = vqmid
  vq = zeros(8,20);
  hump = [12 18 12];
  entry = 1;

  % single hump

  for i=1:18
    vq(entry++,i:i+2) = hump;
  end

  % dual hump

  for i=1:18
    for j=i+3:18
      vq(entry,i:i+2) = hump;
      vq(entry++,j:j+2) = hump;
    end
  end

  figure(1); clf; mesh(vq)

endfunction


function vq = vqmidhigh
  vq = -6*ones(1,30);
  hump = [12 18 12];
  entry = 1;

  % single hump

  for i=1:25
    vq(entry++,i:i+2) = hump;
  end

  % dual hump

  for i=1:25
    for j=i+3:25
      vq(entry,i:i+2) = hump;
      vq(entry++,j:j+2) = hump;
    end
  end

  % triple hump
  
  for i=1:25
    for j=i+5:25
      for k=j+5:25
        vq(entry,i:i+2) = hump;
        vq(entry,j:j+2) = hump;
        vq(entry++,k:k+2) = hump;
      end
    end
  end

  printf("entries: %d\n", entry-1);
  %figure(1); clf; mesh(vq)

endfunction


function vq = vqfull
  K = 20;
  vq = zeros(1,K);
  formant = [12 18 12];
  entry = 1;

  % flat for like background noise and stuff

  vq(entry++,:) = 0;
  
  % single formant, only really likely low down

  for i=1:10
    vq(entry++,i:i+2) = formant;
  end

  % two formants

  for i=1:10
    for j=i+5:K
      vq(entry,i:i+2) = formant;
      vq(entry++,j:j+2) = formant;
    end
  end

  % three formants

  for i=1:10
    for j=i+5:K
      for k=j+5:K
        vq(entry,i:i+2) = formant;
        vq(entry,j:j+2) = formant;
        vq(entry++,k:k+2) = formant;
      end
    end
  end

  printf("entries: %d\n", entry-1);
  %figure(1); clf; mesh(vq)

endfunction

function vq = vqhigh
  vq = zeros(7,10);

  vq(2,:) = [ 0   0   0   0   0  -6 -12 -18 -24 -30];
  vq(3,:) = [ 12 12   6  -6 -12 -18 -24 -30 -36 -42];
  vq(4,:) = [ 6  12  12   6  -6 -12 -18 -24 -30 -36];
  vq(5,:) = [ 0   6  12  12   6  -6 -12 -18 -24 -30];
  vq(6,:) = [ 0   0   6  12  12   6  -6 -12 -18 -24];
  vq(7,:) = [ 0   0   0   6  12  12   6  -6 -12 -18];
  figure(1); clf; plot(vq')
endfunction

