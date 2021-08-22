% trellis.m
% David Rowe July 2021
%
% Testing trellis decoding of Codec 2 Vector Quantiser (VQ)
% information.  Uses soft decision information, probablility of state
% transitions, and left over redundancy to correct errors on VQ
% reception.
%
% VQ indexes are transmitted as codewords mapped to +-1
%
%   y = c + n
%
% where c is the transmitted codeword, y is the received codeword,
% and n is Gaussian noise.
%
% This script generates the test data files:
%
%  cd codec2/build_linux
%  ../script/train_trellis.sh
%

1;

% converts a decimal value to a soft dec binary value
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

    c = -1 + 2*c;
endfunction


% y is vector of received soft decision values (e.g +/-1 + noise)
function [txp indexes] = ln_tx_codeword_prob_given_rx_codeword_y(y, nstates, C)
  nbits = length(y);
  np    = 2.^nbits;
  
  % Find log probability of all possible transmitted codewords
  txp = C * y';
  
  % return most probable codewords (number of states to search)
  [txp indexes] = sort(txp,"descend");  
endfunction

% A matrix of all possible tx codewords C, one per row
function C = precompute_C(nbits)
  np    = 2.^nbits;

  C = zeros(np, nbits);
  for r=0:np-1
    C(r+1,:) = dec2sd(r,nbits);
  end
  
endfunction


% work out transition probability matrix, given lists of current and next
% candidate codewords

function tp = calculate_tp(vq, sd_table, h_table, indexes_current, indexes_next, verbose)
  ntxcw = length(indexes_current);
  tp = zeros(ntxcw, ntxcw);
  for txcw_current=1:ntxcw
    index_current = indexes_current(txcw_current);
    for txcw_next=1:ntxcw
      index_next = indexes_next(txcw_next);
      dist = vq(index_current,:) - vq(index_next,:);
      sd = dist * dist';
      p = prob_from_hist(sd_table, h_table, sd);
      if bitand(verbose, 0x2)
        printf("index_current: %d index_next: %d sd: %f p: %f\n", index_current, index_next, sd, p);
      end
      tp(txcw_current, txcw_next) = log(p);
    end
  end
endfunction


% y is the sequence received soft decision codewords, each row is one
% codeword in time.  sd_table and h_table map SD to
% probability.  Returns the most likely transmitted VQ index c in the
% middle of the codeword sequence y.  We search the most likely ntxcw
% tx codewords out of 2^nbits possibilities.

function c = find_most_likely_index(y, vq, C, sd_table, h_table, nstages, ntxcw, verbose)
    [ncodewords nbits] = size(y);

    % populate the nodes of the trellis with the most likely transmitted codewords

    txp = zeros(nstages, ntxcw); indexes = zeros(nstages, ntxcw);
    for s=1:nstages
      [atxp aindexes] = ln_tx_codeword_prob_given_rx_codeword_y(y(s,:), ntxcw, C);
      txp(s,:) = atxp;
      indexes(s,:) = aindexes;
    end
    
    if verbose
      printf("rx_codewords:\n");
      for r=1:ncodewords
        for c=1:nbits
          printf("%7.2f", y(r,c));
        end
        printf("\n");
      end

      printf("\nProbability of each tx codeword index/binary/ln(prob):\n");
      printf("     ");
      for s=1:nstages
        printf("Time n%+d       ", s - (floor(nstages/2)+1));
      end
      printf("\n");

      for i=1:ntxcw
        printf("%d   ", i);
        for s=1:nstages
	  ind = indexes(s,i) - 1;
          printf("%2d %3s %5.2f   ", ind, dec2bin(ind,nbits), txp(s, i));
        end
        printf("\n");
      end
      printf("\n");
    end
    
    % determine probability of each path
    % (prob of start tx cw)(prob of transition)(prob of next tx cw)
    % cycle through all possibilities of tx cw

    if verbose
      printf("Evaulation of all possible paths:\n");
      printf("  ");
      for s=1:nstages
        printf(" n%+d", s - (floor(nstages/2)+1));
      end
      printf("   indexes");
      printf(" ");
      
      for s=1:nstages
        printf("  txp(%d) ", s-1);
        if s < nstages
          printf(" tp(%d,%d) ", s-1,s);
        end
      end
      printf("     prob  max_prob\n");
    end

    % Determine transition probability matrix for each stage, this
    % changes between stages as lists of candidate tx codewords
    % changes
 
    tp = zeros(nstages, ntxcw, ntxcw);
    for s=1:nstages-1
      tp(s,:,:) = calculate_tp(vq, sd_table, h_table, indexes(s,:), indexes(s+1,:), verbose);
    end

    % OK lets search all possible paths and find most probable

    n = ones(1,nstages); % current node at each stage through trellis, describes current path
    max_prob = -100;
    do
      
      if verbose
       printf("  ");
       for s=1:nstages
          printf("%4d", n(s)-1);
        end
        printf("  ");
         for s=1:nstages
          printf("%2d ", indexes(s,n(s))-1);
        end
      end

      % find the probability of current path
      prob = 0;
      for s=1:nstages
         prob += txp(s, n(s));
         if verbose
           printf("%8.2f ", txp(s, n(s)));
         end
         if s < nstages
	   prob += tp(s, n(s), n(s+1));
           if verbose
             printf("%8.2f ", tp(s, n(s), n(s+1)));
           end
         end
      end

      if (prob > max_prob)
        max_prob = prob;
        max_n = n; 
      end    
    
      if verbose
        printf("%9.2f %9.2f\n", prob, max_prob);
      end

      % next path

      s = nstages;
      n(s)++;
      while (s && (n(s) == (ntxcw+1))) 
        n(s) = 1;
        s--;
        if s > 0
          n(s)++;
        end
      end
    until (sum(n) == nstages)

    c = max_n((nstages+1)/2) - 1;
    if verbose
      printf("\nMost likely path through nodes... ");
      for s=1:nstages
        printf("%4d", max_n(s)-1);
      end
      printf("\nMost likely path through indexes: ");
      for s=1:nstages
        printf("%4d", indexes(s,max_n(s))-1);
      end
      printf("\nMost likely VQ index at time n..: %4d\n", c);
    end
endfunction


% Given a normalised histogram, estimate probability from SD
function p = prob_from_hist(sd_table, h_table, sd)
  p = interp1 (sd_table, h_table, sd, "extrap");
endfunction


% Calculate a normalised histogram of the SD of adjacent frames from
% a file of output vectors from the VQ.
function [sd_table h_table] = vq_hist(vq_output_fn, dec=1)
  K=20; K_st=2+1; K_en=16+1;
  vq_out = load_f32(vq_output_fn, K);
  [r c]= size(vq_out);
  diff = vq_out(dec+1:end,K_st:K_en) - vq_out(1:end-dec,K_st:K_en);
  sd_adj = var(diff');
  [h_table sd_table] = hist(sd_adj,50,1);
endfunction


% vector quantise a sequence of target input vectors, returning the VQ indexes and
% quantised vectors target_
function [indexes target_] = vector_quantiser(vq, target, verbose=1)
  [vq_size K] = size(vq);
  [ntarget tmp] = size(target);
  target_ = zeros(ntarget,K);
  indexes = zeros(1,ntarget);
  for i=1:ntarget
    best_e = 1E32;
    for ind=1:vq_size
      e = sum((vq(ind,:)-target(i,:)).^2);
      if verbose printf("i: %d ind: %d e: %f\n", i, ind, e), end;
      if e < best_e
        best_e = e;
	best_ind = ind;
      end
    end
    if verbose printf("best_e: %f best_ind: %d\n", best_e, best_ind), end;
    target_(i,:) = vq(best_ind,:); indexes(i) = best_ind;
  end
endfunction


% faster version of vector quantiser
function [indexes target_] = vector_quantiser_fast(vq, target, verbose=1)
  [vq_size K] = size(vq);
  [ntarget tmp] = size(target);
  target_ = zeros(ntarget,K);
  indexes = zeros(1,ntarget);

  % pre-compute energy of each VQ vector
  vqsq = zeros(vq_size,1);
  for i=1:vq_size
    vqsq(i) = vq(i,:)*vq(i,:)';
  end

  % use efficient matrix multiplies to search for best match to target
  for i=1:ntarget
    best_e = 1E32;
    e = vqsq - 2*(vq * target(i,:)');
    [best_e best_ind] = min(e);
    if verbose printf("best_e: %f best_ind: %d\n", best_e, best_ind), end;
    target_(i,:) = vq(best_ind,:); indexes(i) = best_ind;
  end
endfunction


% Run trellis decoder over a sequence of frames
function [rx_indexes] = run_test(tx_indexes, vq, sd_table, h_table, ntxcw, nstages, EbNo, verbose)
  [frames nbits]    = size(tx_codewords);
  nerrors           = 0;
  nerrors_simple    = 0;
  tbits             = 0;

  C = precompute_C(nbits);
  % construct tx symbol codewords from VQ indexes
  tx_codewords = zeros(frames, nbits);
  for i=1:frames
    tx_codewords = dec2sd(tx_indexes(i, nbits));
  end  
  rx_codewords = tx_codewords + randn(frames, nbits)/EbNo;
  rx_indexes = zeros(1,frames);
  
  ns2 = floor(nstages/2);
  for f=ns2+1:frames-ns2
    tx_bits        = tx_codewords(f,:) < 0;
    rx_indexes(f)  = find_most_likely_index(rx_codewords(f-ns2:f+ns2,:)/EbNo,
                                            vq, C, sd_table, h_table, nstages, ntxcw, verbose);
    rx_bits        = dec2sd(rx_indexes(f), nbits) < 0;
    rx_bits_simple = rx_codewords(f,:) < 0;
    rx_indexes_van dec_simple_hat(f) = sum(rx_bits_simple .* 2.^(nbits-1:-1:0));

    nerrors        += sum(xor(tx_bits, rx_bits));
    nerrors_simple += sum(xor(tx_bits, rx_bits_simple));
    if verbose
      printf("[%d] %d %d\n", f, nerrors, nerrors_simple);
    end
    tbits += nbits;
  end

  EbNodB = 10*log10(EbNo);
  printf("Eb/No: %3.2f dB nerrors %d %d BER: %3.2f %3.2f std dev: %3.2f %3.2f\n", 
         EbNodB, nerrors, nerrors_simple, nerrors/tbits, nerrors_simple/tbits,
         std(dec_tx_codewords - dec_hat), std(dec_tx_codewords - dec_simple_hat));
  
endfunction

% Simulations ---------------------------------------------------------------------

function test_trellis_against_vanilla(vq_fn, target_fn)
  K = 20; K_st=2+1; K_en=16+1;

  % load VQ
  vq = load_f32(vq_fn, K);
  [vq_size tmp] = size(vq);
  vq = vq(:,K_st:K_en);
  
  % load sequence of target vectors we wish to VQ
  target = load_f32(target_fn, K);

  % lets just test with the first ntarget vectors
  ntarget = 100;
  target = target(1:ntarget,K_st:K_en);
  
  % mean SD of vanilla decode
  [indexes target_ ] = vector_quantiser_fast(vq, target, verbose=0);
  diff = target - target_;
  mse_vanilla = mean(diff(:).^2)

  
endfunction

% Plot histograms of SD at different decimations in time
function vq_hist_dec(vq_output_fn)
  figure(1); clf;
  [sd_table h_table] = vq_hist(vq_output_fn, dec=1);
  plot(sd_table, h_table, "b;dec=1;");
  hold on;
  [sd_table h_table] = vq_hist(vq_output_fn, dec=2);
  plot(sd_table, h_table, "r;dec=2;");
  [sd_table h_table] = vq_hist(vq_output_fn, dec=3);
  plot(sd_table, h_table, "g;dec=3;");  
  [sd_table h_table] = vq_hist(vq_output_fn, dec=4);
  plot(sd_table, h_table, "c;dec=4;");  
  hold off;
  axis([0 300 0 0.5])
  xlabel("SD dB*dB"); title('Histogram of SD(n,n+1)');
endfunction

% Automated tests for vanilla and fast VQ search functions
function test_vq(vq_fn)
  K=20;
  vq = load_f32(vq_fn, K);
  vq_size = 100;
  target = vq(1:vq_size,:);
  indexes = vector_quantiser(target,target, verbose=0);
  assert(indexes == 1:vq_size);
  printf("Vanilla OK!\n");
  indexes = vector_quantiser_fast(target,target, verbose=0);
  assert(indexes == 1:vq_size);
  printf("Fast OK!\n");
endfunction

% Single point test, print out tables
function test_single
  nstages  = 3;
  verbose  = 1;
  nbits = 2;
  ntxcw = 4;
  EbNo = 1;
  
  tx_codewords = [-1 -1; -1 -1; -1 -1];
  rx_codewords = [-1 -1; -1 -1; -1 -1] + randn(nstages, nbits)/EbNo;
  vq = [0 0 0 1;
        0 0 1 0;
	0 1 0 0;
	1 0 0 0];
  sd_table = [0 1 2 4];
  h_table = [0.5 0.25 0.15 0.1];
  C = precompute_C(nbits);
  c = find_most_likely_index(EbNo*rx_codewords, vq, C, sd_table, h_table, nstages, ntxcw, verbose);
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
  nstages  = 3;
  var      = 0.25;
  verbose  = 0;
  frames   = length(codewords_test);
  %frames   = 100;

  tx_codewords = zeros(frames, nbits);
  for f=1:frames
    dec = codewords_test(f, model_param);
    tx_codewords(f,:) = dec2sd(dec, nbits);
  end

  test_traj(tx_codewords, tp(1:nstates,1:nstates,model_param), nstages, var, verbose);

  figure(3)
  plot(codewords_test(1:frames, model_param))

endfunction


% writes an output bitsream file with errors corrected
% replay with something like:
%   $ ./c2dec 700 ../../octave/ve9qrp_10s_trellis.bit - --softdec --natural | play -t raw -r 8000 -s -2 -

function process_test_file(bitstream_filename)
  [bits_per_frame bit_fields] = codec2_700_bit_fields;

  % load training data

  printf("loading training database and generating tp .... ");
  [tp codewords_train] = build_tp("hts.bit");
  printf("done\n");

  % load test data

  printf("loading test database .... ");
  [tp_test codewords_test] = build_tp(bitstream_filename);
  printf("done\n");

  nstages  = 3;
  var      = 0.5;
  verbose  = 0;
  frames   = length(codewords_test);
  %frames   = 150;

  % set up output frames

  frames_out = zeros(1,frames*bits_per_frame);
  frames_out_simple = zeros(1,frames*bits_per_frame);
 
  % run test data through model and trellis decoder

  % just pass these fields straight through

  %for mp=1:10
  for mp=[1 2 3 10]
    nbits = bit_fields(mp);

    % pack parameters in SD format into output frame

    for i=1:frames      
      en = (i-1)*bits_per_frame + sum(bit_fields(1:mp));
      st = en - (nbits-1);
      %printf("fr st: %d offset: %d st: %d en: %d\n", (i-1)*bits_per_frame, sum(bit_fields(1:mp)), st, en);
      frames_out(st:en) = dec2sd(codewords_test(i, mp), nbits);
      frames_out_simple(st:en) = dec2sd(codewords_test(i, mp), nbits);
      %p += bit_fields(i);
    end
  end

  
  % trellis decode these fields

  for mp=4:9
    nbits    = bit_fields(mp);
    nstates  = 2.^nbits;

    printf("processing parameter: %d nbits: %d\n", mp, nbits);

    % normalise transition probability rows

    tpmp = tp(1:nstates,1:nstates, mp);
    s = sum(tpmp+0.001,2);
    [r c] = size(tpmp);
    for i=1:r
      tpmp(i,:) /= s(i);
    end

    tx_codewords = zeros(frames, nbits);
    for f=1:frames
      dec = codewords_test(f, mp);
      tx_codewords(f,:) = dec2sd(dec, nbits);
    end

    [dec_hat dec_simple_hat dec_tx_codewords] = run_test(tx_codewords, tpmp, nstages, var, verbose);

    % pack decoded parameter in SD format into output frame

    for i=1:frames      
      en = (i-1)*bits_per_frame + sum(bit_fields(1:mp));
      st = en - (nbits-1);
      frames_out(st:en) = dec2sd(dec_hat(i), nbits);
      frames_out_simple(st:en) = dec2sd(dec_simple_hat(i), nbits);
    end

  end

  % save frames out

  [tok rem] = strtok(bitstream_filename,".");
  fn_trellis = sprintf("%s_%3.2f_trellis%s",tok,var,rem);
  fn_simple = sprintf("%s_%3.2f_simple%s",tok,var,rem);
  printf("writing output files: %s %s for your listening pleasure!\n", fn_trellis, fn_simple);

  fbitstream = fopen(fn_trellis,"wb");
  fwrite(fbitstream, frames_out,"float32");
  fclose(fbitstream);

  fbitstream = fopen(fn_simple,"wb");
  fwrite(fbitstream, frames_out_simple,"float32");
  fclose(fbitstream);

endfunction


% -------------------------------------------------------------------

graphics_toolkit ("gnuplot");
more off;
randn('state',1);

% uncomment one of the below to run a test or simulation

test_trellis_against_vanilla("../build_linux/vq_stage1.f32", "../build_linux/all_speech_8k_lim.f32")
%test_vq("../build_linux/vq_stage1.f32");
%vq_hist_dec("../build_linux/all_speech_8k_test.f32");
%test_single

