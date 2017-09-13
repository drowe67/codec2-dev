% vq_construct_indep_mag.m
% David Rowe
% Sep 2017
%
%----------------------------------------------------------------------
% Construct a VQ entry using formant prototypes and gain/mag fitting
% indep amplitude and gain for each formant

#{
  ptable - prototype table
         - each row defines a formant prototype with 3 numbers: shape start end
         - shape 1 - gaussian
                 2 - symmetrical triangle
                 3 - triangle with deep fall off at HF side
         - start 0 - one after last prototype center position
                 n - start of search range in 100's of Hz, e.g. 20 -> 2000Hz
         - end   m - end of search range in 100's of Hz
#}

function [idx contrib errors b_log2] = vq_construct_indep_mg(data, w_en = 0)

#{
  ptable = [ 1  4  8;
             2  8 30;
             2 12 30;
             3 30 36];
#}

  ptable = [ 1  4  5 ];
  
  protos = [-4  12 23  24  18  12 8 5 2;
             12 18 12   0   0   0 0 0 0;
             12 18 12 -12 -24 -36 0 0 0];
  protoc = [4 2 2];
  protol = [9 3 6];                      

  nformants = rows(ptable);
  [nRows nCols] = size(data);  
  K = nCols;  
  
  idx = zeros(nRows, nformants); b_log2 = [];

  contrib = zeros(nRows, K);
  
  for r=1:nRows
    % iterate across all formants
    
    % set up target
      
    t = data(r,:);

    % compute weights

    if w_en
      w = find_weights(t, wmin=0.1, tmin=20);
    else
      w = ones(1,K);
    end
    
    figure(10); clf; plot(t); hold on;
    
    for f=1:nformants

      % try formant centre at np points between st and en
      
      shape = ptable(f,1); c_st = ptable(f,2); c_en = ptable(f,3); 
      if c_st == 0  % 0 code means start just after previous formant
        c_st = idx(r,f-1) + 2;
      end
      np = c_en-c_st+1;
      b_log = zeros(np, 2);
      error = zeros(np, 1);

      printf("r: %d f: %d c_st: %d c_en: %d np: %d ------------------------\n", r, f, c_st, c_en, np);

      % test range of centres for formant
      
      for c=c_st:c_en

        % construct current prototype vector and gain window vector

        p_st = c - protoc(shape) + 1; p_en = c - protoc(shape) + protol(shape);
        p = zeros(1,K);
        p(p_st:p_en) = protos(shape, 1:protol(shape));
        g = zeros(1,K);
        g(p_st:p_en) = ones(1, protol(shape));
        
        A = [(p.*w)*p' sum(p.*w); sum(p.*w) sum(g.*w) ];
        d = [(t.*w)*p' (t.*w)*g']';
        b = inv(A)*d;

        v = b(1)*p + b(2)*g;
        diff = t - v;
        plot(v,'g'); plot(diff,'r');
        p = c - c_st + 1;
        b_log(p,:) = b; 
        error(p) = (diff.*w) * diff';

        printf("r: %d f: %d c: %d %d %d e: %f b(1): %f b(2): %f \n", r, f, c, p_st, p_en, error(p), b(1), b(2));
      end

      hold off;
      xx
      % choose best for this formant

      [mn p_min] = min(error);
      % p_min = 1;
      c_min = p_min + c_st - 1;
      b = b_log(p_min,:);
      printf("f: %d cmin: %d b(1): %f b(2): %f -------------\n", f, c_min, b(1), b(2));
      
      % sum this formats contribution
      
      p_st = c_min - protoc(shape) + 1; p_en = c_min - protoc(shape) + protol(shape);
      p = zeros(1,K);
      p(p_st:p_en) = protos(shape, 1:protol(shape));
      g = zeros(1,K);
      g(p_st:p_en) = ones(1, protol(shape));
      v = b(1)*p + b(2)*g;
      contrib(r,:) += v;
      
      % update residual target for next formant search
      
      t = data(r,:) - contrib(r,:);
      
      % log index, mag, gain

      idx(r,f) = c_min;
      b_log2 = [b_log2; b];
    end
    
    errors(r) = sum(t .^ 2); 
  end

endfunction
