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

codec_cmd = "speexenc --bitrate 4000 in.raw - | speexdec - out.raw"

lpc_order = 4;
Fs = 8000;
fo = 200;             % pitch frequency of voice 
wo = 2*pi*fo/Fs;      % pitch in rads
N  = Fs*0.04;         % frame length
frames = Fs/N;        % frames to play through codec
gain = 100;

s = [];
w_log = [];
w__log = [];

for f=1:frames

  % construct LSP symbol

  w = [0.1 0.4  0.5 0.8]*pi;

  % synthesise speech signal

  a=lsptoa(w);
  w_log = [w_log; w];
  l=pi/wo;
  ex = zeros(1,N);
  for m=1:l
    ex += cos(wo*m*(1:N));
  end
  s = [s filter(1, a, ex)];
end

% play through codec

s *= gain;
f=fopen("in.raw","wb");
fwrite(f,s,"short");
fclose(f);
system(codec_cmd);
f=fopen("out.raw","rb");
s_ = fread(f,Inf,"short");
fclose(f);

% extract received symbols from channel and evaluate

for f=1:frames
  a_ = lpcauto(s_((f-1)*N+1:f*N), lpc_order);
  w_ = atolsp(a_);
  w__log = [w__log; w_];
end

figure(1);
clf;
subplot(211)
plot(s)
axis([1 Fs/10 -4000 4000]);
subplot(212)
plot(s_);
axis([1 Fs/10 -4000 4000]);

figure(2)
plot(w_log, 'g', w__log, 'r+')
