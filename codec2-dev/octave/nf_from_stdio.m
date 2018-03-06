% nf_from_gr.m
% David Rowe Mar 2018

#{

  Calculate NF in real time from 16 bit real samples from stdin

  1/ When used with UDP feature of gqrx in USB mode:

    $ nc -ul 7355 | octave --no-gui -qf nf_from_stdio.m

  2/ Using command line tools:

  a) Airspy:
  
    $ airspy_rx -a 6000000 -l 14 -m 15 -v 15 -r - -f 434.998 -b 1 | \
      csdr convert_s16_f | csdr fir_decimate_cc 50 | csdr convert_f_s16 | \
      octave --no-gui -qf ~/codec2-dev/octave/nf_from_stdio.m 120000 complex

  a) HackRF:
  
    Term 1:

    $ ~/codec2-dev/octave$ nc -ul 7355 | octave --no-gui -qf nf_from_stdio.m 80000 complex

    Term 2:
    
    $ hackrf_transfer -r - -f 434995000 -s 4000000 -a 1 -p 1 -l 40 -g 32 | \
      csdr convert_s8_f | csdr fir_decimate_cc 50 | csdr convert_f_s16 | \
      nc localhost -u 735
   
  TODO:
    [ ] work out why noise power st bounces around so much, signal power seems stable
    [ ] reduce CPU load, in particular of plotting
#}

graphics_toolkit ("gnuplot")

% command line arguments

arg_list = argv ();
if nargin == 0
  printf("\nusage: %s FsHz [real|complex] [testToneLeveldBm]\n\n", program_name());
  exit(0);
end

Fs = str2num(arg_list{1});
shorts_per_sample = 1;

if nargin == 2
  if strcmp(arg_list{2}, "real")
    shorts_per_sample = 1;
  end
  if strcmp(arg_list{2}, "complex")
    shorts_per_sample = 2;
  end
end
  
Pin_dB = -100; % level of input test tone
if nargin == 3
  Pin_dB = str2num(arg_list{3});
end

printf("Fs: %d shorts_per_sample: %d Pin_dB: %f\n", Fs, shorts_per_sample, Pin_dB);

[s,c] = fread(stdin, shorts_per_sample*Fs, "short");

while c
  if shorts_per_sample == 2
    s = s(1:2:end)+j*s(2:2:end);
  end
  S = fft(s.*hanning(Fs));
  SdB = 20*log10(abs(S));
  figure(1); plot(real(s)); axis([0 Fs -3E4 3E4]);
  figure(2); plot(SdB); axis([0 12000 40 160]);

  % assume sine wave is between 2000 and 4000 Hz, and dominates energy in that
  % region.  Noise is between 5000 - 10000 Hz
  
  sig_st = 2000; sig_en = 5000;
  noise_st = 6000; noise_en = 10000;

  % find peak and sum power a few bins either side, this ensure we don't capture
  % too much noise as well
  
  [pk pk_pos] = max(abs(S));
  if pk_pos > 5
    Pout_dB1 = 10*log10(sum(abs(S(pk_pos-5:pk_pos+5)).^2));     % Rx output power with test signal
  else
    Pout_dB1 = 0;
  end
  
  Pout_dB = 10*log10(sum(abs(S(sig_st:sig_en)).^2));      % Rx output power with test signal
  G_dB    = Pout_dB - Pin_dB;                             % Gain of Rx                           
  Nout_dB = 10*log10(sum(abs(S(noise_st:noise_en)).^2)/(noise_en-noise_st));  % Rx output power with noise 
  Nin_dB  = Nout_dB - G_dB;                               % Rx input power with noise
  No_dB   = Nin_dB; %- 10*log10(noise_en-noise_st);       % Rx input power with noise in 1Hz bandwidth
  NF_dB   = No_dB + 174;                                  % compare to thermal noise to get NF
  printf("Pout: %4.1f  %d %4.1f Nout: %4.1f G: %4.1f  No: %4.1f NF: %3.1f dB\n", Pout_dB, pk_pos, Pout_dB1, Nout_dB, G_dB, No_dB, NF_dB);

  pause(2);
  [s,c] = fread(stdin, shorts_per_sample*Fs, "short");
endwhile

