% LDPC unit test script
% David Rowe 18 Dec 2013
% Based on siulation by Bill Cowley

% Start CML library

currentdir = pwd;
addpath '/home/david/tmp/cml/mat'    % assume the source files stored here
cd /home/david/tmp/cml
CmlStartup                            % note that this is not in the cml path!
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

code_param = ldpc_init(rate, framesize, modulation, mod_order, mapping);

Ntrials = 84;
EsNo=10;

Tbits = Terrs = Ferrs = 0;
    
data = [];
r = []; 

% Encode a bunch of frames

for nn = 1: Ntrials        
    d = round( rand( 1, code_param.data_bits_per_frame ) );
    data = [data d];
    [codeword, s] = ldpc_enc(d, code_param);
    code_param.code_bits_per_frame = length(codeword);
    code_param.symbols_per_frame = length(s);
    r = [r s];
end

% Decode a bunch of frames

for nn = 1: Ntrials        
    st = (nn-1)*code_param.symbols_per_frame + 1;
    en = (nn)*code_param.symbols_per_frame;
    detected_data = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, r(st:en), EsNo);
    st = (nn-1)*code_param.data_bits_per_frame + 1;
    en = (nn)*code_param.data_bits_per_frame;
    error_positions = xor( detected_data(1:code_param.data_bits_per_frame), data(st:en) );
    Nerrs = sum( error_positions);
        
    if Nerrs>0, fprintf(1,'x'),  else fprintf(1,'.'),  end
    if (rem(nn, 50)==0),  fprintf(1,'\n'),  end    
    if Nerrs>0,  Ferrs = Ferrs +1;  end
    Terrs = Terrs + Nerrs;
    Tbits = Tbits + code_param.data_bits_per_frame;        
end
fprintf(1,'\n')

