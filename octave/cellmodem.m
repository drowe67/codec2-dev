% cellmodem.m
% David Rowe May 2005
%
% Simulation of modem for low rate date through a cell phone codec
%
% Ideas:
%   + insert low rate codec
%   + generate bunch of symbols, run through codec, measure MSE, choose best set
%   + measure probablility of error, "distance" from other symbols
%   + can we use VQ training algorithm for this?  Start with 2 symbols,
%     pass through channel, measure MSE, split again?
%   + start with cmd line version of codec, frame synchronous
%   + add symbol timing estimator later
%   + try different frame rates
%   + simulate impairments like HP/LP filtering.  Can we correct for this?
%   + pilots symbols so we can use energy as well?
%   + set of F0 as well
%   + LSP quantisers preferrentially preserve peaks, so three peaks is a
%     reasonable signal set.  Modulate position and bandwidth to creat 
%     symbol set.  Maybe 50 or 100 Hz grid for LSPs, evaluate that some how.

graphics_toolkit ("gnuplot");
lsp;

lpc_order = 4;
Fs = 8000;
fo = 100;             % pitch frequency of voice 
wo = 2*pi*fo/Fs;      % pitch in rads
N  = Fs*0.04;         % frame length

% construct LSP symbol

w = [0.1 0.4  0.5 0.8]*pi

% synthesise speech signal

a=lsptoa(w)
l=pi/wo;
ex = zeros(1,N);
for m=1:l
    ex += cos(wo*m*(1:N));
end
s = filter(1, a, ex);

% insert codec here (e.g. write to file, run codec, read file)

% extract received symbol from channel

a_ = lpcauto(s, lpc_order)
w_ = atolsp(a_)

figure(1);
clf;
subplot(211)
plot(s)
subplot(212)
plot(w,'g+',w_','r+');
axis([0 lpc_order+1 0 pi])
