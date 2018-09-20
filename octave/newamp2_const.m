% newamp2_const.m
% Constants for newamp2 candidate C prototype codec

Fs      = 8000;   % sample rate
max_amp = 160;    % maximum number of rate L samples (time varying basedon Wo)
K       = 20;     % number of samples after interpolation to fixed rate K
M       = 2;      % decimation rate from codec processing frame rate

Nbitspercodecframe          = 56;
Nprotectedbitspercodecframe = 32;
Nunprotect                  = Nbitspercodecframe - Nprotectedbitspercodecframe;

