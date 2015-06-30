% trellis.m
% David Rowe June 2015
%
% Testing ideas for trellis decoding of codec parameters.  Uses
% soft decision information, probablility of state transitions,
% and left over redundancy to correct errors on codec parameter
% reception.
%
% [ ]  Add probability of transition
% [ ]  Will this mess up a perfectly received signal, i.e. all
%      codewords perfectly received.  This should dominate prob
%      calc (with small No)

graphics_toolkit ("gnuplot");
more off;

function [tp codewords] = build_tp(bitstream_filename)
  bits_per_frame = 28; 
  bit_fields = [1 5 3 3 2 4 3 3 2 2];
  nfields = length(bit_fields);

  % load encoded speech file (one float per bit)

  fbitstream = fopen(bitstream_filename, "rb"); 
  bitstream = fread(fbitstream, "float32");
  fclose(fbitstream);
  nframes = floor(length(bitstream)/bits_per_frame);
  bitstream = bitstream < 0;
  
  % extract LSPs codewords

  codewords = zeros(nframes,nfields);
  for f=1:nframes
    aframe = bitstream((f-1)*bits_per_frame+1:f*bits_per_frame);
    field = 1;
    st = bit_fields(field);
    for l=1:nfields
      nbits = bit_fields(field);
      %printf("st: %d %d\n", st, st+nbits-1);
      codeword = aframe(st:st+nbits-1)';
      %printf("codeword: %d\n", codeword);
      codewords(f,l) = sum(codeword .* 2.^(nbits-1:-1:0));
      %printf("nbits: %d lsp(%d, %d) = %d\n", nbits, f, l, codewords(f,l));
      st += nbits;
      field++;
    end
  end

  % determine transition probablilities of each codeword
  % state at time n, state a time n+1
  % prob must add to 1
  % so for every state, of every codeword, compute histogram of 
  % transition probabilities

  tp = zeros(32,32,nfields);
  for f=1:nframes-1
    for l=1:nfields
      acodeword = codewords(f,l);
      bcodeword = codewords(f+1,l);
      %printf("%d %d %d %d\n", f, l, acodeword, bcodeword);
      tp(acodeword+1,bcodeword+1,l)++;
    end
  end
endfunction


% y is vector of +/- 1 soft decision values for 0,1 transmitted bits

function lnp = ln_prob_of_tx_codeword_c_given_rx_codeword_y(y)
  nbits = length(y);
  np    = 2.^nbits;

  % work through all possible received codeworks and determine probability
  % given a number of bits

  lnp = zeros(1,np);
  for i=0:np-1

    % convert to binary

    c = zeros(1,nbits);
    for j=0:nbits-1
      mask = 2.^j;
      if bitand(i,mask)
        c(j+1) = 1;
      end
    end

    % map to +/- 1

    c = 1 - 2*c;
    
    % probability calculation for this i

    lnp(i+1) = sum(y .* c);
  end

endfunction


% y is the received soft decision codedwords, each row is one codeword in time
% tp is the transition probabilities, each row is the start state
% returns the most likely transitted codeword c

function c = find_most_likely_codeword(y, tp)
    [nstates nbits] = size(y);
    tp
    y_n   = y(1,:);  % codeword received at time n
    y_n_1 = y(2,:);  % codeword received at time n+1

    lnp_n   = ln_prob_of_tx_codeword_c_given_rx_codeword_y(y_n);
    lnp_n_1 = ln_prob_of_tx_codeword_c_given_rx_codeword_y(y_n_1);

    % determine probability of a path
    % (prob of start state)(prob of transition)(prob of next state)
    % cycle through all possibilities of states
    % state variable for possible state
    % step through them exhaustively

    max_prob = 0;
    for n=1:nstates
      for n_1=1:nstates
        prob = lnp_n(n);
        prob += tp(n, n_1);
        prob += lnp_n_1(n_1);
        printf("state: %d %d %f\n", n, n_1, prob);
        if (prob > max_prob)
          max_prob = prob;
          max_n = n; max_n_1 = n;
        end
      end
    end
    
endfunction

% for a given parameter
% weight for No
% given received symbols y(n) and y(n+1)
% calculate all possible paths
% select most likely


printf("loading training database and generating tp .... ");
[tp codewords] = build_tp("hts.bit");
printf("done\n");

y = [-1 -1; -1 -1];
find_most_likely_codeword(y,tp(1:4,1:4,5))

% TODO: normalise tp (is this rqd? Just a const?), ln(tp), add No to calculation
%       then test!  Try to insert in channel and correct errors
