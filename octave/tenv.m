% tenv.m
% David Rowe August 2022 & Feb 2023
%
% Exploring time domain envelope & phase algorithms to minimise PAPR

1;

function fg = do_plots(fg, h, H, F0, Wo, L, fftx)
  figure(fg++); clf;
  subplot(211); plot(fftx, 20*log10(abs(h)));
  hold on; plot((1:L)*F0,20*log10(abs(H)),'r+','markersize',20); hold off;
  subplot(212); plot(fftx, angle(h));
  hold on; plot((1:L)*F0,angle(H),'r+','markersize',20); hold off;  
 
  % synthesise speech

  N = 320;
  s = zeros(1,N);
  for m=1:L
    s = s + abs(H(m))*cos(m*Wo*(1:N) + angle(H(m)));
  end
  figure(fg++); clf;
  papr_dB = 10*log10(max(abs(s).^2)/mean(abs(s).^2));
  l = sprintf("b-;CPAPR: %3.1f dB;", papr_dB);
  plot(s,l);
endfunction

Fs = 8000;

% set up second order system
w = 2*pi*600/Fs; gamma = 0.0;
ak = [1 -2*gamma*cos(w) gamma*gamma];

% sample magnitude and phase using a large number of points
Nfft = 512;
h = freqz(1,ak,Nfft/2+1);
fftx = (0:Nfft/2)*(Fs/Nfft);

% sample at harmonics of F0
F0=100; Wo = 2*pi*F0/Fs; L = floor(Fs/(2*F0));
H = freqz(1,ak,(1:L)*Wo);

fg = 1;
fg = do_plots(fg, h, H, F0, Wo, L, fftx);

% use iterative method to determine phase spectra to minimise PAPR

theta = zeros(1,L);
A = abs(H);
P = Fs/F0;

for m=2:L
  % synthesise one pitch period with phases up to m-1
  s = zeros(1,P);
  for k=1:m-1
     s = s + A(k)*cos(k*Wo*(0:P-1) + theta(k));
  end
   % determine current peak position
  [max_s max_n] = max(abs(s));
  % set phase of m-th to minimise peak
  if sign(max(s(max_n))) > 0
    theta(m) = -pi-m*Wo*(max_n-1);
  else
    theta(m) = -m*Wo*(max_n-1);
  end
  figure(5); clf; hold on; plot(s); plot(A(m)*cos(m*Wo*(0:P-1) + theta(m))); hold off;
  printf("m: %d max_n: %d sign(s(max_n)) %f theta(m): %f\n", m, max_n, sign(s(max_n)), theta(m));
  fflush(stdout);
  %k = kbhit();
end

H1 = A .* exp(j*theta);
fg = do_plots(fg, h, H1, F0, Wo, L, fftx);

 
