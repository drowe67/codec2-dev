% hackrf_uc.m
%
% David Rowe Nov 2015
%
% Upconverts a real baseband sample file to a file suitable for input into a HackRF
%
% To play file at 10.7MHz used:
%   octave:25> hackrf_uc("fsk_10M.iq","fsk_horus_rx_1200_96k.raw")
%   $ hackrf_transfer -t ../octave/fsk_10M.iq -f 10000000 -a 0 -a 1 -x 40

function hackrf_uc(outfilename, infilename)
  Fs1 = 96E3;  % input sample rate
  Fs2 = 10E6;  % output sample rate to HackRF
  fc = 700E3;  % offset to shift to, HackRF doesn't like signals in the centre
  A  = 100;    % amplitude of signal after upc-nversion (max 127)
  N  = Fs1*20;
  
  fin = fopen(infilename,"rb");
  s1 = fread(fin,Inf,"short");
  fclose(fin);
  ls1 = length(s1);

  % limit noise to first 4 kHz.  Otherwise noise dominates signals and
  % gain control below pushes wanted signal amplitude down so far we
  % use only a few DAC/ADC bits

  [b a] = cheby2(6,40,[500 4000]/(Fs1/2));
  %s1 = filter(b,a,s1);

  % single sided freq shifts, we don't want DSB

  s1 = hilbert(s1(1:N)); 

  % upsample to Fs2

  M = Fs2/Fs1;
  s2 = resample(s1(1:N),Fs2,Fs1);
  ls2 = length(s2);
  mx = max(abs(s2));
  t = 0:ls2-1;

  % shift up to Fc, not sure of rot90 rather than trasnpose operator '
  % as we don't want complex conj, that would shift down in freq

  sout = rot90((A/mx)*s2) .* exp(j*2*pi*t*fc/Fs2);

  save_hackrf(outfilename,sout);
end
