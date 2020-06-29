% Test FSK at symbol level in CML using RA LDPC code,  Bill, June 2020
% Initialise CML first; 
% If required setup numb of codewords (Ncw)   and channel bits/sym (bps)
% and use plt to control debug plots 

ldpc;

%rand('seed',1);
%randn('seed',1);
format short g 
more off
init_cml('~/cml/');



Hfilename = 'HRAa_1536_512.mat'
load(Hfilename);
[Nr Nc] = size(HRA);
Nbits = Nc - Nr;
Krate = (Nc-Nr)/Nc
framesize = Nc;
[H_rows, H_cols] = Mat2Hrows(HRA);
code_param.H_rows = H_rows;
code_param.H_cols = H_cols;
code_param.P_matrix = [];

if exist('Ncw')==0,  Ncw=100,   end
if exist('plt')==0,    plt=0;   end
if exist('bps')==0,    bps=4;   end 


M=2^bps;   nos =0;    clear res
S =CreateConstellation('FSK', M);
if M==2,   Ebvec=[7:0.2: 8.5],   end 
if M==4,   Ebvec=[4:0.3: 7],  end 
if M==16,  Ebvec=[2.5:0.25:4.4];   end 


disp(['Symbol-based ' num2str(M) 'FSK sim with rate 3/4 code, ' ...
       num2str(Ncw) ' 2k codewords'])

for Eb = Ebvec
    
    Ec =  Eb + 10*log10(Krate);
    Es =  Ec + 10*log10(bps);
    Eslin = 10^(Es/10);       %Es/N0 = 1/2k_n^2
    
    Terrs =0;
    for nn = 1:Ncw
        
        txbits = randi(2,1,Nbits) -1;
        
        codeword = LdpcEncode( txbits, code_param.H_rows, code_param.P_matrix );
        code_param.code_bits_per_frame = length( codeword );
        code_param.data_bits_per_frame = length(txbits);
        Nsymb = code_param.code_bits_per_frame/bps;      
        
        Tx = Modulate(codeword, S);
        
        kn = sqrt(1/(2*Eslin));
        Rx = Tx + kn * (randn(size(Tx)) + j*randn(size(Tx)));
        
        SNRlin = Eslin;  % Valenti calls this SNR, but seems to be Es/N0
        symL = DemodFSK(Rx, SNRlin, 2);    %demod type is nonCOH, without estimate amplitudes
        bitL = Somap(symL);
        if plt>0, figure(110);   hist(bitL);   title('bit LLRs')
            figure(111);   hist(bitL);   title('Sym Ls'),  pause, 
        end
        max_it =100;   decoder_type =0;
        
        [x_hat, PCcnt] = MpDecode( -bitL, code_param.H_rows, code_param.H_cols, ...
            max_it, decoder_type, 1, 1);
        Niters = sum(PCcnt~=0); 
        detected_data = x_hat(Niters,:);
        error_positions = xor( detected_data(1:code_param.data_bits_per_frame), txbits );
        
        Nerrs = sum( error_positions);
        if plt>1, figure(121);  plot(error_positions);  Nerrs,  end 
        Terrs = Terrs + Nerrs; 
    end
    
    BER = Terrs/ (Ncw*Nbits);
    
    %HDs = (sign(bitL)+1)/2;
    %NerrsHD  = sum(codeword~=HDs);
    %BER_HD = Nerrs/Nbits;
    
    nos = nos+1;
    res(nos, :) = [Eb, Terrs,  BER]
    
end
figure(91)
semilogy(res(:,1), res(:,3)); grid on; hold on;
%semilogy(res(:,1), berawgn(res(:,1), 'fsk', M, 'noncoherent'), 'g');
title('MFSK BER with LDPC FEC')





