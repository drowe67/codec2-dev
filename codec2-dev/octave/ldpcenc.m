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
Eprob = 0.15;

vocoderframesize = 52;
nvocoderframes = 8;
nbitspermodemframe = 72;

code_param = ldpc_init(rate, framesize, modulation, mod_order, mapping);

data = [];
r = []; 

% Encoder: Generate simulated vocoder data, insert UW, and LPDC encode ---------------

Nframes = 100;
uw = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];

% repeat same simulated vocoder data to ease testing

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

printf("Encoded %d LDPC frames\n", Nframes);

% Modulator: Modulate to QPSK symbols ------------------------------------------

lpackedcodeword=length(packedcodeword);
fc=fopen("codeword.bin","rb");
fm=fopen("modcodeword.bin","wb");
lpackedmodem = nbitspermodemframe/8;
n = 0;

[packedmodem, count] = fread(fc,lpackedmodem,"uchar");
while (count == lpackedmodem)
    n++;
    unpackedmodem = unpackmsb(packedmodem);

    ii = 1;
    for i=1:2:length(unpackedmodem)
        mod_unpackedmodem(ii) = qpsk_mod(unpackedmodem(i:i+1));
        mod_unpackedmodem_float32(i) = real(mod_unpackedmodem(ii));
        mod_unpackedmodem_float32(i+1) = imag(mod_unpackedmodem(ii));
        ii += 1;
    end

    fwrite(fm, mod_unpackedmodem_float32, "float32");
    [packedmodem, count] = fread(fc,lpackedmodem,"uchar");
end
fclose(fc);
fclose(fm);
printf("Modulated %d modem frames\n", n);


% Decoder: Sync with LDPC frames, LDPC decode, strip off UW, measure BER -------

fm=fopen("modcodeword.bin","rb");

mod_uw = build_mod_uw(uw, 2*length(vd)/length(uw));

mod_codeword = zeros(1, code_param.code_bits_per_frame/2);
lmod_codeword = code_param.code_bits_per_frame/2;

Terrs = 0; Ferrs = 0; Tbits = 0; Tframes = 0; nerr = [];
corr = []; n = 0;
sync_state = 0; sync_count = 0;

[mod_unpackedmodem_float32, count] = fread(fm,nbitspermodemframe, "float32");
while (count == nbitspermodemframe)
    n++;

    mod_unpackedmodem = mod_unpackedmodem_float32(1:2:nbitspermodemframe) + j*mod_unpackedmodem_float32(2:2:nbitspermodemframe);
    erasures = rand(1,length(mod_unpackedmodem)) < Eprob; 
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

    % read in one modulated modem frame at a time

    [mod_unpackedmodem_float32, count] = fread(fm, nbitspermodemframe, "float32");
end

fprintf(1,"\nFrames: %d bits: %d errors: %d BER = %f FER = %f\n", Tframes, Tbits, Terrs, Terrs/Tbits, Ferrs/Tframes);
subplot(211)
plot(corr);
subplot(212)
plot(nerr);
