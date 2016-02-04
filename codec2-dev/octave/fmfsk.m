%
% fmfsk.m
% Author: Brady O'Brien 3 Feb 2016
%   Copyright 2016 David Rowe
%  
%  All rights reserved.
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License veRbion 2, as
%  published by the Free Software Foundation.  This program is
%  distributed in the hope that it will be useful, but WITHOUT ANY
%  WARRANTY; without even the implied warranty of MERCHANTABILITY or
%  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
%  License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License
%  along with this program; if not, see <http://www.gnu.org/licenses/>.

% mancyfsk.m modem, extracted and made suitable for C implementation


1;

% Init fmfsk modem
%Fs is sample frequency
%Rb is pre-manchester bit rate
function states = fmfsk_init(Fs,Rb)
    assert(floor(Fs/(Rb*2))==(Fs/(Rb*2)));
    assert(mod(Fs,Rb*2)==0);
    
    states.Rb = Rb;
    states.Rs = Rb*2;   % Manchester-encoded bitrate
    states.Fs = Fs;
    states.Ts = Fs/states.Rs;
    % Processing buffer size. about 40ms here
    states.N = floor(states.Rs*.040)*states.Ts;     
    states.nin = states.N;          % Samples in the next demod cycle
    states.nstash = states.Ts*2;    % How many samples to stash away between proc cycles for timing adjust
    states.nmem =  states.N+(4*states.Ts);
    states.nsym = floor(states.Rs*.040)
    %Old sample memory
    %states.oldsamps = zeros(1,states.nstash);
    
    states.oldsamps = zeros(1,states.nmem);
    
    %Last sampled-stream output, for odd bitstream generation
    states.lastint = 0;
    
    %Some stats
    states.norm_rx_timing = 0;
    
endfunction

%Generate a stream of manchester-coded bits to be sent
% to any ordinary FM modulator or VCO or something
function tx = fmfsk_mod(states,inbits)
    Ts = states.Ts;
    tx = zeros(1,length(inbits)*2);
    for ii = 1:length(inbits)
        st = 1 + (ii-1)*Ts*2;
        md = st+Ts-1;
        en = md+Ts;
        if inbits(ii)==0
            tx(st:md)   = -ones(1,Ts);
            tx(md+1:en) =  ones(1,Ts);
        else
            tx(st:md)   =  ones(1,Ts);
            tx(md+1:en) = -ones(1,Ts);
        end
    end
endfunction

%Demodulate a bag of bits from the output of an FM demodulator
% This function produces nbits output bits and takes states.nin samples
function [rx_bits states] = fmfsk_demod(states,rx)
    Ts = states.Ts;
    Fs = states.Fs;
    Rs = states.Rs;
    nin = states.nin;
    N = states.N;
    nsym = states.nsym;
    nbits = states.nsym/2;
    nmem = states.nmem;
    nstash = states.nstash;
    
    nold = nmem-nin;
    ssamps = states.oldsamps;
    
    %Shift in nin samples
    ssamps(1:nold) = ssamps(nmem-nold+1:nmem);
    ssamps(nold+1:nmem) = rx;
    stats.oldsamps = ssamps;
    
    rx_filt = zeros(1,nsym*Ts);
    %Integrate Ts input samples at every offset
    %This is the same thing as filtering with a filter of all ones
    % out to Ts.
    % It's implemented like this for ease of C-porting
    for ii=(1:(nsym+1)*Ts)
        st = ii;
        en = st+Ts-1;
        rx_filt(ii) = sum(ssamps(st:en));
    end
 
    % Fine timing estimation ------------------------------------------------------

    % Estimate fine timing using line at Rs/2 that Manchester encoding provides
    % We need this to sync up to Manchester codewords. 
    Np = length(rx_filt);
    w = 2*pi*(Rs)/Fs;
    x = (rx_filt .^ 2) * exp(-j*w*(0:Np-1))';
    norm_rx_timing = angle(x)/(2*pi) - 0.42
    rx_timing = round(norm_rx_timing*Ts)
    
    %If rx timing is too far out, ask for more or less sample the next time
    % around to even it all out
    next_nin = N;
    if norm_rx_timing > 0.25
       next_nin += Ts/2;
    end
    if norm_rx_timing < -0.25;
       next_nin -= Ts/2;
    end
    states.nin = next_nin;
    
    % Sample rx_filt at the optimum inst, as figured by rx_timing
    rx_symsamp = rx_filt((Ts/2)+Ts+rx_timing:Ts:length(rx_filt));
    length(rx_symsamp)
    
    rx_bits = zeros(1,nbits);
      
      
endfunction

