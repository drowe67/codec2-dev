function llrs = mfsk_hd_to_llrs(sd, map, SNR=4)
% compare these results to SD method of mfsk_sd_to_llrs
% Bill June 2020


  [x, ind] = max(sd');     % find biggest SD
  b2 = map(ind', :)';      % get bit pairs
  b1 = b2(:)';             % covert to row vector
  b1 = b1*2 -1;            % convert to -1/ +1
  llrs = -4*SNR*b1;         % rough scaling;  default is ~6dB SNR
  
endfunction
