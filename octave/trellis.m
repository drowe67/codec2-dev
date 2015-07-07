% trellis.m
% David Rowe June 2015
%
% Testing ideas for trellis decoding of codec parameters.  Uses
% soft decision information, probablility of state transitions,
% and left over redundancy to correct errors on codec parameter
% reception.
%
% [X]  Add probability of transition
% [X]  Will this mess up a perfectly received signal, i.e. all
%      codewords perfectly received.  This should dominate prob
%      calc (with small No)
% [ ] can we measure erroneous decodes?  Whn it makes it worse?
%     + like all FEC, at some point, it will get a poorer BER
%     + can we drawa curve of coded versus uncoded errors?
% [ ] test with actual data, add errors, see if it corrects any
% [ ] can we draw a trellis with actual values?
%
% Note we need SD bits and natural binary order encoding, for example
%   ./c2enc 700 ../../raw/ve9qrp_10s.raw ../../octave/ve9qrp_10s.bit --softdec --natural

graphics_toolkit ("gnuplot");
more off;
randn('state',1);

% Before I couldn't even sp3ll functional programmer.  Now I are one ;)

function [bits_per_frame bit_fields] = codec2_700_bit_fields
  bits_per_frame = 28;                % number of bits/frame for "700" mode
  bit_fields = [1 5 3 3 2 4 3 3 2 2]; % number of bits in each field for "700" mode
  % voiced, Wo, energy, LSP 1..6, 2 spare
endfunction

% builds a metric of trasnition probablilities for each codec field (model parameter)

function [tp codewords] = build_tp(bitstream_filename)
  [bits_per_frame bit_fields] = codec2_700_bit_fields;
  nfields = length(bit_fields);

  % load encoded speech file (one float per bit)

  fbitstream = fopen(bitstream_filename, "rb"); 
  bitstream = fread(fbitstream, "float32");
  fclose(fbitstream);
  nframes = floor(length(bitstream)/bits_per_frame);
  bitstream = bitstream < 0;
  
  % extract each field (model parameter) for bit stream

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


% converts a decimal vlaue to a soft dec binary value

function c = dec2sd(dec, nbits)
    
    % convert to binary

    c = zeros(1,nbits);
    for j=0:nbits-1
      mask = 2.^j;
      if bitand(dec,mask)
        c(nbits-j) = 1;
      end
    end

    % map to +/- 1

    c = 1 - 2*c;
endfunction

% y is vector of +/- 1 soft decision values for 0,1 transmitted bits

function lnp = ln_prob_of_tx_codeword_c_given_rx_codeword_y(y)
  nbits = length(y);
  np    = 2.^nbits;

  % work through all possible received codeworks and determine probability
  % given a number of bits

  lnp = zeros(1,np);
  for i=0:np-1

    c = dec2sd(i,nbits);

    % probability calculation for this i

    lnp(i+1) = sum(y .* c);
  end

endfunction


% y is the received soft decision codedwords, each row is one codeword in time
% tp is the transition probabilities, each row is the start state
% returns the most likely transitted codeword c

function c = find_most_likely_codeword(y, tp, verbose)
    [ncodewords nbits] = size(y);
    nstates = 2^nbits;
    y_n   = y(1,:);  % codeword received at time n
    y_n_1 = y(2,:);  % codeword received at time n+1
    y_n_2 = y(3,:);  % codeword received at time n+1

    lnp_n   = ln_prob_of_tx_codeword_c_given_rx_codeword_y(y_n);
    lnp_n_1 = ln_prob_of_tx_codeword_c_given_rx_codeword_y(y_n_1);
    lnp_n_2 = ln_prob_of_tx_codeword_c_given_rx_codeword_y(y_n_2);

    if verbose
      % printf received codewords

      printf("rx_codewords\n");
      for r=1:ncodewords
        for c=1:nbits
          printf("%7.2f", y(r,c));
        end
        printf("\n");
      end

      % print probability of each rx codeword table

      printf("\nProbability of each tx codeword:\n");
      printf("tx_codeword   lnp(n)  lnp(n+1)  lnp(n+2)  \n");
      printf("bin dec\n");
      for i=1:nstates
        printf("%3s  %2d      %7.2f   %7.2f   %7.2f\n", dec2bin(i-1,nbits), i-1, lnp_n(i), lnp_n_1(i), lnp_n_2(i));
      end
      printf("\n");

      printf("\nTransition Probability table:\n");
      printf("state/next_state\n");
      for r=1:nstates
        for c=1:nstates
          printf("%7.2f", tp(r,c));
        end
        printf("\n");
      end
      printf("\n");
    end

    % determine probability of each path
    % (prob of start state)(prob of transition)(prob of next state)
    % cycle through all possibilities of states
    % state variable for possible state
    % step through them exhaustively

    if verbose
      printf("Evaulation of all possible paths\n");
      printf("tx codeword\n");
      printf(" n  n+1  n+2    lnp(n) lnp(n+1) lnp(n+2)  tp(n,n+1)  tp(n+1,n+2)      total\n");
    end
    max_prob = -100;
    for n=1:nstates
      for n_1=1:nstates
        for n_2=1:nstates
          prob = lnp_n(n);
          prob += lnp_n_1(n_1);
          prob += lnp_n_2(n_2);
          prob += tp(n, n_1);
          prob += tp(n_1, n_2);
          if verbose
            printf("%2d %4d %4d  %7.2f  %7.2f  %7.2f     %6.2f       %6.2f   %8.2f\n", n-1, n_1-1, n_2-1, lnp_n(n), lnp_n_1(n_1), lnp_n_2(n_2), tp(n, n_1), tp(n_1, n_2), prob);
          end
          if (prob > max_prob)
            max_prob = prob;
            max_n = n; max_n_1 = n_1; max_n_2 = n_2;
          end 
        end
      end
    end
    
    if verbose
      printf("And the winner is.....\n");
      printf("%2d %4d %4d  %7.2f  %7.2f  %7.2f     %6.2f       %6.2f   %8.2f\n", max_n-1, max_n_1-1, max_n_2-1, lnp_n(n), lnp_n_1(max_n_1), lnp_n_2(max_n_2), tp(max_n, max_n_1), tp(max_n_1, max_n_2), max_prob);
    end

    c = max_n_1 - 1;

endfunction

% Simulations ---------------------------

function [dec_hat dec_simple_hat] = run_test(tx_codewords, transition_probabilities, var, verbose)
  [frames nbits]    = size(tx_codewords);
  nerrors           = 0;
  nerrors_simple    = 0;
  tbits             = 0;
  nstates           = 2.^nbits;

  rx_codewords = tx_codewords + sqrt(var)*randn(frames, nbits);
  dec_hat = dec_simple_hat = zeros(1,frames);
  
  for f=2:frames-1
    tx_bits        = tx_codewords(f,:) < 0;
    dec_hat(f)     = find_most_likely_codeword(rx_codewords(f-1:f+1,:)/(2*var),log(transition_probabilities+0.001), verbose);
    rx_bits        = dec2sd(dec_hat(f), nbits) < 0;
    rx_bits_simple = rx_codewords(f,:) < 0;
    dec_simple_hat(f) = sum(rx_bits_simple .* 2.^(nbits-1:-1:0));

    nerrors        += sum(xor(tx_bits, rx_bits));
    nerrors_simple += sum(xor(tx_bits, rx_bits_simple));
    if verbose
      printf("[%d] %d %d\n", f, nerrors, nerrors_simple);
    end
    tbits += nbits;
  end

  EbNodB = 10*log10(1/(2*var));
  printf("Eb/No: %3.2f dB nerrors %d %d BER: %3.2f %3.2f\n", 
         EbNodB, nerrors, nerrors_simple, nerrors/tbits, nerrors_simple/tbits);
  
endfunction


% Single point test, print out tables

function test_single
  tp      = [1 0 0 0; 1 0 0 0; 1 0 0 0; 1 0 0 0];
  var     = 0.001;
  verbose = 1;

  tx_codewords = [-1 1; 1 -1; -1 1];
  run_test(tx_codewords, tp, var, verbose)
endfunction


% plot a trajectory of codewords, to help test "distance" is better using trellis than simple decoding

function test_traj(tx_codewords, tp, var, verbose)
  [frames nbits] = size(tx_codewords);
  nstates = 2.^nbits;

  tp = tp ./ sum(tp+0.001,2);

  [dec_hat dec_simple_hat] = run_test(tx_codewords, tp, var, verbose);

  dec_tx_codewords = zeros(1,frames);
  for f=2:frames-1
    tx_bits = tx_codewords(f,:) < 0;
    dec_tx_codewords(f) = sum(tx_bits .* 2.^(nbits-1:-1:0));
  end

  figure(1);
  subplot(211)
  stem(dec_tx_codewords - dec_hat)
  axis([1 frames -nstates/2 nstates/2]);
  title('Trellis Decoding');
  subplot(212)
  stem(dec_tx_codewords - dec_simple_hat);
  axis([1 frames -nstates/2 nstates/2]);
  title('Simple Decoding');

  figure(2)
  tp
  length(1:nstates)
  mesh(1:nstates,1:nstates,tp)
  ylabel('current state');
  xlabel('next state');
endfunction


% A contrived trajectories to check out the idea

function simple_traj
  %tp      = [1 0 0 0; 1 0 0 0; 1 0 0 0; 1 0 0 0];
  tp      = [0.5 0.25 0 0; 0.25 0.5 0.25 0; 0 0.25 0.5 0.25; 0 0 0.25 0.5];
  var     = 0.5;
  verbose = 0;
  nbits   = 2;
  frames  = 100;

  tx_codewords = ones(frames, nbits);
  test_traj(tx_codewords, tp, var, verbose)
endfunction


% Run through a trajectory, say from a different file to what we
% trained with.  Plot trajectory before and after passing through
% our decoder.  Try first with no noise, then with noise.

function test_codec_model_parameter(bitstream_filename, model_param)
  [bits_per_frame bit_fields] = codec2_700_bit_fields;

  % load training data

  printf("loading training database and generating tp .... ");
  [tp codewords_train] = build_tp("hts.bit");
  printf("done\n");

  % load test data

  printf("loading test database .... ");
  [tp_test codewords_test] = build_tp(bitstream_filename);
  printf("done\n");

  % run test data through

  nbits    = bit_fields(model_param);
  nstates  = 2.^nbits;
  var      = 0.5;
  verbose  = 0;
  frames   = length(codewords_test);
  frames   = 100;

  tx_codewords = zeros(frames, nbits);
  for f=1:frames
    dec = codewords_test(f, model_param);
    tx_codewords(f,:) = dec2sd(dec, nbits);
  end

  test_traj(tx_codewords, tp(1:nstates,1:nstates,model_param), var, verbose);

  figure(3)
  plot(codewords_test(1:frames, model_param))

  figure(4);
  fs=fopen("../raw/ve9qrp_10s.raw","rb");
  s = fread(fs,Inf,"short");
  plot(s(1:320*frames))
  axis([1 320*frames -3E4 3E4 ])
endfunction


% uncomment one of the below to run a simulation

%test_single
%simple_traj;
test_codec_model_parameter("ve9qrp_10s.bit", 6);


