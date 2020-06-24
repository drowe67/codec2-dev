% Plot some results from "fsk_llr_test.m"
% After running, array "res" contains rows of simulation results  
%   Eb  Ec  M  Ltype  Nbits  Nerr BERraw   
%   (some uncoded rows contains -1 to indicate val is not applicable)
 
figure(102);   hold on;    

%uncoded results
sub = res(res(:,4)==-1 & res(:,3)==2, :) 
semilogy(sub(:,1), sub(:,7), 'k+--')
sub = res(res(:,4)==-1 & res(:,3)==4, :) 
semilogy(sub(:,1), sub(:,7), 'k--')

% coded results 
for M = [2 4 ] 
lt = ''; 
if M==2
   lt = '-+'
   sub = res(res(:,4)==1 & res(:,3)==M, :) 
   semilogy(sub(:,1), sub(:,6)./sub(:,5), ['k' lt])
endif 

sub = res(res(:,4)==2 & res(:,3)==M, :)
semilogy(sub(:,1), sub(:,6)./sub(:,5), ['g' lt])

sub = res(res(:,4)==3 & res(:,3)==M, :)
semilogy(sub(:,1), sub(:,6)./sub(:,5), ['b' lt])
endfor

ylabel('BER')
xlabel('Eb/N0 (Info Bits; dB)')
title('MFSK LLR test (solid is coded; +is 2FSK, grn is SDprob, blk is SDorig, blu is HD)')
print   -dpng LLRsMFSK.png

