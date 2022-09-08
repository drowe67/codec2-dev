% rateK_batch.m
%
% David Rowe 2022
%
% batch mode rateK resampling
%
% Usage:
%   Make sure codec2-dev is compiled with the -DDUMP option - see README.md for
%    instructions.
%   ~/codec2-dev/build_linux/src$ ./c2sim ../../raw/hts1a.raw --dump hts1a
%   $ cd ~/codec2-dev/octave
%   octave:14> rateK_batch("../build_linux/src/hts1a")

function rateK_batch(samname)
  more off;
  
  [Espline F0] = arateK_batch(samname,'spline');
  Epara = arateK_batch(samname,'para');
  figure(1); clf;
  plot(F0,Espline,'b+;spline;'); hold on;
  plot(F0,Epara,'g+;para;'); xlabel('F0 (Hz)'); hold off;
  xlabel('F0 (Hz)'); ylabel('E dB^2');
  figure(2); clf;
  [hspline nnspline] = hist(Espline,50);
  [hpara nnpara] = hist(Epara,50); hold on;
  title('Histogram of E');
endfunction

function [E F0] = arateK_batch(samname, resampler='spline')
  more off;
  
  newamp_700c;
  Fs = 8000; K = 80;

  model_name = strcat(samname,"_model.txt");
  model = load(model_name);
  [frames tmp] = size(model);

  [rate_K_surface sample_freqs_kHz] = resample_const_rate_f_lin(model, K, resampler);
  model_ = resample_rate_L(model, rate_K_surface, sample_freqs_kHz, resampler);
  fn = strcat(samname,"_out_model.txt");
  save("-ascii",fn,"model_");

  E = zeros(1,frames);
  for f=1:frames
    Wo = model(f,1); F0(f) = Fs*Wo/(2*pi); L = model(f,2);
    Am  = model_(f,3:(L+2)); AmdB  = 20*log10(Am);
    Am_ = model(f,3:(L+2)); AmdB_ = 20*log10(Am_);
    Lmin = round(200/F0(f)); Lmax = floor(3700/F0(f));
    E(f)  = sum((AmdB(Lmin:Lmax) - AmdB_(Lmin:Lmax)).^2)/(Lmax-Lmin+1);
  end
  printf("%8s mean SD: %4.2f dB^2\n", resampler, mean(E));
endfunction

