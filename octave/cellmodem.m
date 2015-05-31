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
rand('state',1);
lsp;

% given a vector of LSP symbols, constructs a synthesised speech signals
% anfd runs it through a codec

function [w__log mse] = run_sim(sim_in, w_log)
  N         = sim_in.N;
  Wo        = sim_in.Wo;
  frames    = sim_in.frames;
  lpc_order = sim_in.lpc_order;

  s      = [];
  w__log = [];
  L      = floor(pi/Wo);
  phi    = zeros(1,L);

  for f=1:frames

    % synthesise speech signal

    a=lsptoa(w_log(f,:));

    ex = zeros(1,N);
    for m=1:L
      phi(m) += Wo*m*N;
      ex += cos(phi(m) + Wo*m*(0:N-1));
    end

    s = [s filter(1, a, ex)];
  end

  % play through codec

  s *= sim_in.gain;
  f=fopen("in.raw","wb");
  fwrite(f,s,"short");
  fclose(f);
  system(sim_in.codec_cmd);
  f=fopen("out.raw","rb");
  s_ = fread(f,Inf,"short");
  fclose(f);

  % extract received symbols from channel and evaluate

  mse = zeros(1,frames);
  for f=1:frames
    a_ = lpcauto(s_((f-1)*sim_in.N+1:f*N), lpc_order);
    w_ = atolsp(a_);
    w__log = [w__log; w_];
    error = w__log(f,:) - w_log(f,:);
    mse(f) = error*error';
  end
endfunction

% constants -----------------------------------------------------

sim_in.codec_cmd = "speexenc --bitrate 4000 in.raw - | speexdec - out.raw";

lpc_order   = sim_in.lpc_order = 6;
Fs          = sim_in.Fs = 8000;
              sim_in.Fo = 100;             % pitch frequency of voice 
              sim_in.Wo = 2*pi*Fo/Fs;      % pitch in rads
N           = sim_in.N  = Fs*0.04;         % frame length
frames      = sim_in.frames = 1000;        % frames to play through codec
              sim_in.gain = 100;
Nsym        = 8;                           % number of symbols to find

% start with some LSP random vectors
% for stable filter most be monotonically increasing on 0..pi

w_log = [];
for f=1:frames;
  w = sort(rand(1,lpc_order)*pi);
  w_log = [w_log; w];
end
[w__log mse] = run_sim(sim_in, w_log);

% sort by MSE to get the best symbols

[sort_mse sort_ind] = sort(mse);
symbols = w_log(sort_ind(1:Nsym),:)

% Play these symbols through the codec in random order

w_log = [];
symb_ind = [];
for f=1:frames
  symb_ind(f) = floor(1 + rand(1,1)*Nsym);
  w_log = [w_log; symbols(symb_ind,:)];
end

[w__log mse] = run_sim(sim_in, w_log);

% now see if we can "detect" them

for f=1:frames

  % check received symbol against codebook of symbols

  min_e = 1E6;
  for i=1:Nsym
    e = w__log(f,:)*symbols(i,:)';
    if e < min_e
      min_e = e;
      min_ind = i;
    end
  end

  symb_ind_out(f) = min_ind;
end

figure(1)
clf
plot(symb_ind,'g',symb_ind_out,'r');
