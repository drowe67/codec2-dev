% ldpcenc.m
% David Rowe 20 Dec 2013
% 
% LDPC encoder test program. Encodes and modulates a random data stream 

% Start CML library

currentdir = pwd;
addpath '/home/david/tmp/cml/mat'    % assume the source files stored here
cd /home/david/tmp/cml
CmlStartup                           % note that this is not in the cml path!
cd(currentdir)

% Our LDPC library

ldpc;

% Start simulation

rate = 3/4; 
framesize = 576;  

mod_order = 4; 
modulation = 'QPSK';
mapping = 'gray';

demod_type = 0;
decoder_type = 0;
max_iterations = 100;

vocoderframesize = 52;
nvocoderframes = 8;

code_param = ldpc_init(rate, framesize, modulation, mod_order, mapping);

data = [];
r = []; 

% Encode a bunch of frames

Nframes = 100;

% repeat same codeword frame for now to ease testing

vd = round( rand( 1, vocoderframesize*nvocoderframes) );
d = insert_uw(vd, [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]);

data = [data d];
[codeword, s] = ldpc_enc(d, code_param);
code_param.code_bits_per_frame = length(codeword);
code_param.symbols_per_frame = length(s);
packedcodeword = packmsb(codeword);

fc=fopen("codeword.bin","wb");
for nn = 1: Nframes        
    fwrite(fc,packedcodeword,"char");
end
fclose(fc);

printf("framesize: %d data_bits_per_frame: %d code_bits_per_frame: %d\n", ...
        framesize, code_param.data_bits_per_frame,  code_param.code_bits_per_frame);


 
