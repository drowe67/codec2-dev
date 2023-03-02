% tenv.m
% David Rowe August 2022 & Feb 2023
%
% Exploring time domain envelope & phase algorithms to minimise PAPR

1;

% use iterative method to determine phase spectra to minimise PAPR
function theta = phase_from_mag_iterative(A,Fs,F0,verbose=0)
  Wo = 2*pi*F0/Fs; L = floor(Fs/(2*F0)); P = Fs/F0;

  theta = zeros(1,L);
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
    if verbose
     printf("m: %d max_n: %d sign(s(max_n)) %f theta(m): %f\n", m, max_n, sign(s(max_n)), theta(m));
    end
  end
  
  % map to -pi .. pi
  theta -= 2*pi*round(theta/(2*pi));
  
endfunction

function s = synthesise_speech(A,theta,Wo)
  L = floor(pi/Wo);
 
  N = 320;
  s = zeros(1,N);
  for m=1:L
    s = s + A(m)*cos(m*Wo*(1:N) + theta(m));
  end
endfunction

% run simulation using a second order system
function fg = experiment1(fg, title_str, gamma)
  Fs = 8000;
  
  % set up second order system
  w = 2*pi*600/Fs;
  ak = [1 -2*gamma*cos(w) gamma*gamma];

  % sample magnitude and phase using a large number of points
  Nfft = 512;
  h = freqz(1,ak,Nfft/2+1);
  fftx = (0:Nfft/2)*(Fs/Nfft);
  
  % sample at harmonics of F0
  F0=100; Wo = 2*pi*F0/Fs; L = floor(Fs/(2*F0));
  H = freqz(1,ak,(1:L)*Wo);
  A = abs(H); theta = angle(H);
  theta_iter = phase_from_mag_iterative(abs(H),Fs,F0);
  
  figure(fg++); clf;
  subplot(211); plot(fftx, 20*log10(abs(h)),'b'); title(title_str);
  hold on; plot((1:L)*F0,20*log10(A),'b+','markersize',20); hold off;
  subplot(212); plot(fftx, angle(h));
  hold on; plot((1:L)*F0,theta_iter,'r+','markersize',20); hold off; 
   
  s = synthesise_speech(A,theta,Wo);
  s_iter = synthesise_speech(A,theta_iter,Wo);
  papr_dB = 10*log10(max(abs(s).^2)/mean(abs(s).^2));
  l = sprintf("b-;CPAPR: %3.1f dB;", papr_dB);
  figure(fg++); clf;
  subplot(211); plot(s,l); title(title_str);
  papr_dB = 10*log10(max(abs(s_iter).^2)/mean(abs(s).^2));
  l = sprintf("r-;iter CPAPR: %3.1f dB;", papr_dB);
  subplot(212); plot(s_iter,l);
endfunction

fg = 1;
fg = experiment1(fg,'Voiced',0.95);
fg = experiment1(fg,'Unvoiced',0);

