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
EsNo = 10;

vocoderframesize = 52;
nvocoderframes = 8;

code_param = ldpc_init(rate, framesize, modulation, mod_order, mapping);

data = [];
r = []; 

% Encode a bunch of frames

Nframes = 100;
uw = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];

% repeat same codeword frame for now to ease testing

vd = round( rand( 1, vocoderframesize*nvocoderframes) );
d  = insert_uw(vd, uw);

data = [data d];
[codeword, s] = ldpc_enc(d, code_param);
code_param.code_bits_per_frame = length(codeword);
code_param.symbols_per_frame = length(s);
packedcodeword = packmsb(codeword);

fc=fopen("codeword.bin","wb");
for nn = 1: Nframes        
    fwrite(fc,packedcodeword,"uchar");
end
fclose(fc);

%printf("framesize: %d data_bits_per_frame: %d code_bits_per_frame: %d\n", ...
%        framesize, code_param.data_bits_per_frame,  code_param.code_bits_per_frame);

% rx simulation (separate later)

mod_uw = build_mod_uw(uw, 2*length(vd)/length(uw));

lpackedcodeword=length(packedcodeword);
fc=fopen("codeword.bin","rb");
lpackedmodem = 72/8;
mod_codeword = zeros(1, code_param.code_bits_per_frame/2);
lmod_codeword = code_param.code_bits_per_frame/2;

for m=1:8

    % read in one modem frame at a time

    packedmodem = fread(fc,lpackedmodem,"uchar");
    unpackedmodem = unpackmsb(packedmodem);

    j = 1;
    for i=1:2:length(unpackedmodem)
        mod_unpackedmodem(j) = qpsk_mod(unpackedmodem(i:i+1));
        j += 1;
    end

    % keep buffer of one entire codeword

    mod_codeword(1:lmod_codeword-length(mod_unpackedmodem)) = mod_codeword(length(mod_unpackedmodem)+1:lmod_codeword);
    mod_codeword(lmod_codeword-length(mod_unpackedmodem)+1:lmod_codeword) = mod_unpackedmodem;

    uw_sync = look_for_uw(mod_codeword(1:length(mod_uw)), mod_uw);
    if (uw_sync)
        % force UW symbols as they are known (is this needed?)

        % LDPC decode

        detected_data = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, mod_codeword, EsNo);

        % unpack payload data, removing UW

        vd_rx = remove_uw(detected_data(1:code_param.data_bits_per_frame), length(vd), length(uw));

        % measure BER

        error_positions = xor(vd, vd_rx);
        Nerrs = sum(error_positions);
        if Nerrs>0, fprintf(1,'x'),  else fprintf(1,'.'),  end

        % save packed payload data to disk
    end
end

fprintf(1,'\n')
fclose(fc);

