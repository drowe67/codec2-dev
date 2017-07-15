% c2wideband_const.m
% David Rowe July 2017
%
% Constants for wideband Codec 2, bit like #define for C

Fs = 16000;  % sampel rate in Hz
K = 30;      % number of Mel spaced amplitude samples

Tf = 0.01;   % frame period in seconds

dec = 2;     % decimation factor.  10ms update of the core
             % sinusuodial codec is rather fast and reducing it to
             % 20ms makes very little difference in quality, but
             % halves the data going into the DCT/quantisation.
             % While the DCT should encode highly correlated data
             % like 10ms frame efficiently, there is some noise in
             % the core paraneters estimation so nice to use a lower
             % frame rate if possible.

Nt = 8;      % number of blocks in time.  Trade off between latency
             % and coding efficiency.  If the DCT has larger blocks to
             % play with it can remove more correlation, up to a limit
             % where the data is no longer correlated in time or freq.
             % Also set to match LDPC code frame size

