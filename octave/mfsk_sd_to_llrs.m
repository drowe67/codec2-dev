function llrs = mfsk_sd_to_llrs(SD, map, v, sig)
  % Array of MFSK SD is converted to bit LLR vector according to mapping 'map' 
  % SD is M elements wide;  map is M by log2(M); v and sig are params of Rician pdf 
  % Bill (VK5DSP) for Codec2, version 24/6/20 
  
  Smap = size(map); 
  Ssd = size(SD); 
  if Smap(1) == 2
    sig = sig *1.2;    %       <<<   scaling factor TBC 
    llrs = zeros(1, Ssd(1)); 
    bps = 1; 
    assert(2^bps == Ssd(2), 'wrong SD length?')
    % note that when v=0, the pdf reverts to Rayleigh 
    for kk = 1: Ssd(1)
      Prx_1 =  rice(SD(kk,2), v, sig) * rice(SD(kk,1),0,sig);
      Prx_0 =  rice(SD(kk,2), 0, sig) * rice(SD(kk,1),v,sig); 
      Llim = 1e-5; 
      if or(isnan(Prx_1),  Prx_1<Llim), Prx_1 = Llim; end  
      if or(isnan(Prx_0),  Prx_0<Llim), Prx_0 = Llim; end  
      llrs(kk) = -log(Prx_1/Prx_0);
      %% note inversion of LLRs 
   endfor 
  else
    if map==[0 0; 0 1; 1 0; 1 1]   % a specific case for 4FSK
       llrs = zeros(2, Ssd(1));  
       for kk = 1: Ssd(1)
          % pn_m is prob of sub car n given bit m is transmitted 
          p1_0 = rice(SD(kk,1), 0, sig);  p1_1 = rice(SD(kk,1), v, sig);
          p2_0 = rice(SD(kk,2), 0, sig);  p2_1 = rice(SD(kk,2), v, sig);
          p3_0 = rice(SD(kk,3), 0, sig);  p3_1 = rice(SD(kk,3), v, sig);
          p4_0 = rice(SD(kk,4), 0, sig);  p4_1 = rice(SD(kk,4), v, sig);

          %   second bit 
          Prx_1 =  p1_0*p3_0*(p2_1*p4_0 + p2_0*p4_1);
          Prx_0 =  p2_0*p4_0*(p1_1*p3_0 + p1_0*p3_1);

          Llim = 1e-5; 
          if or(isnan(Prx_1),  Prx_1<Llim), Prx_1 = Llim; end  
          if or(isnan(Prx_0),  Prx_0<Llim), Prx_0 = Llim; end  
          llrs(2,kk) = -log(Prx_1/Prx_0);
       
	  %  1st  bit   (LHS) 
          Prx_1 =  p1_0*p2_0*(p3_1*p4_0 + p3_0*p4_1);
          Prx_0 =  p3_0*p4_0*(p1_1*p2_0 + p1_0*p2_1);

          Llim = 1e-5; 
          if or(isnan(Prx_1),  Prx_1<Llim), Prx_1 = Llim; end  
          if or(isnan(Prx_0),  Prx_0<Llim), Prx_0 = Llim; end  
          llrs(1,kk) = -log(Prx_1/Prx_0);
       endfor
       llrs= llrs(:);   % convert to single vector of LLRs
       llrs = llrs'; 
    else
      disp('the mapping is '); map  
      error('not sure how that mapping works!!')
      endif
  endif
endfunction

