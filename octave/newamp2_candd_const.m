% newamp2_candd_const.m
% Constants for newamp2 candidate D prototype codec

Fs      = 8000;   % sample rate
max_amp = 160;    % maximum number of rate L samples (time varying basedon Wo)
K       = 40;     % number of samples after interpolation to fixed rate K
dec     = 2;      % decimation rate from codec processing frame rate

Nbitspercodecframe          = 56;

