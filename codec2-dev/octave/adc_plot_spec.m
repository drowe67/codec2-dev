#!/usr/bin/env octave
% Plot the spectrum and waveform coming off of the ADC

% Author: Brady O'Brien 28 June 2016

%   Copyright 2016 Brady O'Brien
%  
%  All rights reserved.
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License version 2.1, as
%  published by the Free Software Foundation.  This program is
%  distributed in the hope that it will be useful, but WITHOUT ANY
%  WARRANTY; without even the implied warranty of MERCHANTABILITY or
%  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
%  License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License
%  along with this program; if not, see <http://www.gnu.org/licenses/>.


%adc_file_name = '/home/baobrien/workspace/freetel-code/codec2-dev/stm32/adc_samp'
adc_file_name = '/dev/ttyACM0'
%adc_file_name = 'infifo';
fs = 96000
sampsize = 19200
%graphics_toolkit('gnuplot')

fin = fopen(adc_file_name,'r')
first = 1
cont = 0
sampinc = (1/sampsize)*fs

[samps cont] = fread(fin,sampsize,'short');
y = (1:(sampsize/2));
yp = (1:sampsize/2)*(fs/(sampsize/2));
pltdat = plot(10*log10(abs(fft(samps)(y))));
pltx = 10*log10(abs(fft(samps)(y)));
axis([0 fs 30 80])
set (pltdat, "ydatasource", "pltx"); 
set (pltdat, "xdatasource", "yp"); 
%set (pltdat, "ydatasource", "y"); 

while(first || cont==sampsize)
	first = 0;
	%sleep(.001);
	[samps cont] = fread(fin,sampsize,'short');
	pltx = 10*log10(abs(fft(samps)(y)));
	refreshdata
	refresh
endwhile
