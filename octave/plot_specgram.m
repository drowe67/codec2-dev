% plot_specgram.m
% David Rowe May 2017
%
% As the name suggests.....


function S = plot_specgram(x)

  % set step size so we end up with an aray of pixels that is square.  Otherwise
  % imagesc will only use a ribbon on the figure.

  Fs = 8000;
  l = length(x);

  Nfft = 512; Nfft2 = Nfft/2; window = 512; nstep = window/2;
  nsamples = floor(l/nstep) - 1;
  S = zeros(nsamples,Nfft);
  h = hanning(Nfft);
  for i=1:nsamples
    st = (i-1)*nstep+1; en = st+window-1;
    %printf("i: %d st: %d en: %d l: %d\n", i, st, en, l);
    S(i,:) = 20*log10(abs(fft(x(st:en).*h')));
  end
  mx = ceil(max(max(S))/10)*10;
  S = max(S,mx-30);
  mesh((1:Nfft2)*4000/Nfft2, (1:nsamples)*nstep, S(:,1:Nfft2));
  view(90,-90);
  %xlabel('Frequency (Hz)');
  %ylabel('Sample');
en
