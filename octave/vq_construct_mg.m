%----------------------------------------------------------------------
% Construct a VQ entry using formant prototypes and gain/mag fitting

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

function [idx contrib errors b_log2] = vq_construct_mg(data)

  ptable = [ 1  4  8 1 10;
             2  8 30 0 30;
             2 12 30 0 30;
             3 30 36 0 40];
  
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

  #{
    range of F1, specify range of each formant
    specify shape

    [X] table of shapes
    [X] work out centre position
    [X] another table to refer to each shape based on formant number
    [X] const for number of formants
    [ ] iterate over positions
    [ ] gain/mag match at each position, MSE over entire vector
    [ ] way to graphically peer inside what's going on
    [ ] how to model flat in first 1kHz?
        + cld perhaps exclude search entirely as an option
        + or perhap sit will trun up low mags
  #}
  
  for r=1:nRows
    % iterate across all formants
    
    for f=1:nformants

      % try formant centre at np points between st and en
      
      shape = ptable(f,1); c_st = ptable(f,2); c_en = ptable(f,3); 
      if c_st == 0  % 0 code means start just after previous formant
        c_st = idx(r,f-1) + 2;
      end
      np = c_en-c_st+1;
      b_log = zeros(np, 2);
      error = zeros(np, 1);

      % we are matching a target segment between st and en
       
      st = ptable(f,4); en = ptable(f,5);
      if st == 0
        st = c_st;
      end
      lt = en-st+1;
      t = data(r,st:en);

      printf("r: %d f: %d c_st: %d c_en: %d np: %d st: %d en: %d ------------------------\n",
              r, f, c_st, c_en, np, st, en);

      % test range of centres for formant
      
      for c=c_st:c_en
        v = zeros(1, K);
        v_st = c - protoc(shape) + 1; v_en = c - protoc(shape) + protol(shape);
        v(v_st:v_en) = protos(shape, 1:protol(shape));
        v = v(st:en);
        A = [v*v' sum(v); sum(v)  lt];
        d = [t*v' sum(t)]';
        b = inv(A)*d;
        b(1) = max(0.5,b(1));
        diff = t - (b(1)*v + b(2));
        p = c - c_st + 1;
        b_log(p,:) = b; 
        error(p) = diff * diff';

        %printf("r: %d f: %d c: %d e: %f b(1): %f b(2): %f \n", r, f, c, error(p), b(1), b(2));
      end

      % choose best for this formant

      [mn p_min] = min(error);
      c_min = p_min + c_st - 1
      b = b_log(p_min,:);
      v = zeros(1, K);
      v_st = c_min - protoc(shape) + 1; v_en = c_min - protoc(shape) + protol(shape);
      v(v_st:v_en) = protos(shape, 1:protol(shape));
      contrib(r,st:en) = b(1)*v(st:en) + b(2);
      
      % log index, mag, gain, might need a matric with one row per formant

      idx(r,f) = c_min;
      b_log2 = [b_log2; b];
    end
    
    errors(r) = sum((data(r,:) - contrib(r,:)) .^ 2); 
  end

endfunction
