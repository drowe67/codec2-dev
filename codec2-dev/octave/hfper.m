% hfper.m
% David Rowe 2 June 2014
% Quick and dirty HF PER calculator/simulator

function hfper(ber, nbits, ntrials)

  % Raw PER with no FEC --------------------------------------

  nper = 0;
  for i=1:ntrials
      nerr = sum(rand(1,nbits) < ber);
      if nerr >0
          nper++;
      end
  end
  printf("Raw PER..................: %0.3f\n", nper/ntrials);

  % Half rate block code, e.g. Golay (23,12) with 3 bit error
  % correcting capability

  % Golay (23,12) that can correct 3 errors (fails at 4) ------

  ncodeword = 23;
  ncorrect = 3;
  nper = 0;
  for i=1:ntrials
      nerr = sum(rand(1,ncodeword) < ber);
      if nerr > ncorrect
          nper++;
      end
  end
  printf("One Golay codeword.......: %0.3f\n", nper/ntrials);

  % Several Golay codewords concatenated ----------------------

  m = floor(nbits/12);        % number of codewords

  nper = 0;
  for i=1:ntrials

    % test each codeword in packet, if any of the codewords has > 4
    % errors, entire packet is a dud

    no_errors = 1;
    for k=1:m
      nerr = sum(rand(1,ncodeword) < ber);
      if (nerr > ncorrect) && no_errors
          nper++;
          no_errors = 0;
      end
    end

  end
  printf("Packet protected by Golay: %0.3f\n", nper/ntrials);
  
endfunction
 
