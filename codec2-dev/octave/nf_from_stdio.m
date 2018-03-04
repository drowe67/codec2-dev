% nf_from_gr.m
% David Rowe Mar 2018

#{
  Calculate NF in real time from 16 bit real samples from stdin
#}

N = 48000;
Pin_dB = -100;

graphics_toolkit ("gnuplot")

[s,c] = fread(stdin, N, "short");

while c
  S = fft(s.*hanning(N));
  SdB = 20*log10(abs(S));
  figure(1); plot(real(s)); axis([0 N -3E4 3E4]);
  figure(2); plot(SdB); axis([0 12000 40 160]);

  % assume sine wave is between 2000 and 4000 Hz, and dominates energy in that
  % region.  Noise is between 5000 - 10000 Hz
  
  sig_st = 2000; sig_en = 4000;
  noise_st = 5000; noise_en = 10000;
  
  Pout_dB = 10*log10(var(S(sig_st:sig_en)));      % Rx output power with test signal
  G_dB    = Pout_dB - Pin_dB;                     % Gain of Rx                           
  Nout_dB = 10*log10(var(S(noise_st:noise_en)));  % Rx output power with noise
  Nin_dB  = Nout_dB - G_dB;                       % Rx input power with noise
  No_dB   = Nin_dB - 10*log10(noise_en-noise_st); % Rx input power with noise in 1Hz bandwidth
  NF_dB   = No_dB + 174;                          % compare to thermal noise to get NF
  printf("Pout: %4.1f  Nout: %4.1f G: %4.1f  No: %4.1f NF: %3.1f dB\n", Pout_dB, Nout_dB, G_dB, No_dB, NF_dB);

  pause(2);
  [s,c] = fread(stdin, N, "short");
endwhile

