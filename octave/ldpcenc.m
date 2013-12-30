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

rand('state',1);

rate = 3/4; 
framesize = 576;  

mod_order = 4; 
modulation = 'QPSK';
mapping = 'gray';

demod_type = 0;
decoder_type = 0;
max_iterations = 100;
EsNo = 10;
Eprob = 0.18;

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

Terrs = 0; Ferrs = 0; Tbits = 0; Tframes = 0; nerr = [];
corr = []; n = 0;
sync_state = 0; sync_count = 0;

[packedmodem, count] = fread(fc,lpackedmodem,"uchar");
while (count == lpackedmodem)
    n++;
    unpackedmodem = unpackmsb(packedmodem);

    j = 1;
    for i=1:2:length(unpackedmodem)
        mod_unpackedmodem(j) = qpsk_mod(unpackedmodem(i:i+1));
        j += 1;
    end

    erasures = rand(1,length(mod_unpackedmodem))<Eprob; 
    mod_unpackedmodem(erasures) = 0;

    % keep buffer of one entire codeword

    mod_codeword(1:lmod_codeword-length(mod_unpackedmodem)) = mod_codeword(length(mod_unpackedmodem)+1:lmod_codeword);
    mod_codeword(lmod_codeword-length(mod_unpackedmodem)+1:lmod_codeword) = mod_unpackedmodem;

    [uw_sync corr(n)] = look_for_uw(mod_codeword(1:length(mod_uw)), mod_uw);
    if (uw_sync)
      sync_state = 1;
    end

    if (sync_state && (sync_count == 0))
        Tframes++;

        % force UW symbols as they are known (is this needed?)

        % LDPC decode

        detected_data = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, mod_codeword, EsNo);

        % unpack payload data, removing UW

        vd_rx = remove_uw(detected_data(1:code_param.data_bits_per_frame), length(vd), length(uw));

        % measure BER

        error_positions = xor(vd, vd_rx);
        Nerrs = sum(error_positions);
        if Nerrs>0, fprintf(1,'x'); Ferrs++; ,  else fprintf(1,'.'),  end
        Tbits += length(vd);
        Terrs += Nerrs;
        nerr(Tframes) = Nerrs;

        % save packed payload data to disk
    end

    if (sync_state)
        sync_count++;
        if (sync_count == 8)
            sync_count = 0;
        end
    end

    % read in one modem frame at a time

    [packedmodem, count] = fread(fc, lpackedmodem, "uchar");
end
fclose(fc);

fprintf(1,"\nFrames: %d bits: %d errors: %d BER = %f FER = %f\n", Tframes, Tbits, Terrs, Terrs/Tbits, Ferrs/Tframes);
subplot(211)
plot(corr);
subplot(212)
plot(nerr);
