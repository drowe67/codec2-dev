#!/usr/bin/octave -qf

% fsk_horus_stream.m
% David Rowe 13 Oct 2015
%
% Experimental near space balloon FSK demodulator, takes 8kHz 16 bit samples from 
% stdin, output txt string on stdout
%
% usage:
%  $ chmod 777 fsk_horus_stream.m
%  $ rec -t raw -r 8000 -s -2 -c 1 - -q | ./fsk_horus_stream.m
%
% or to test with a stored file:
%  $ cat ~/Desktop/vk5arg-3.wav | ./fsk_horus_stream.m
%

fsk_horus;  % include library (make sure calls to functions at bottom are commented out)

more off;
states = fsk_horus_init();
N = states.N;
Rs = states.Rs;
nsym = states.nsym;
nin = states.nin;
nfield = states.nfield;
npad = states.npad;
uw = states.uw;
EbNo = 0;
SNR = 0;

rx = [];
rx_bits_buf = [];

[s,c] = fread(stdin, N, "short");

while c

  rx = [rx s'];
 
  % demodulate to bit stream

  while length(rx) > nin
    [rx_bits states] = fsk_horus_demod(states, rx(1:nin)');
    rx_bits_buf = [rx_bits_buf rx_bits];
    rx = rx(nin+1:length(rx));
    nin = states.nin;
    EbNo = 0.9*EbNo + 0.1*states.EbNodB;
    SNR = EbNo + 10*log10(states.Rs/3000);
    %printf("nin: %d length(rx): %d length(rx_bits_buf): %d \n", nin, length(rx), length(rx_bits_buf));
  endwhile
  f = (states.f1+states.f2)/2; shift = states.f2 - states.f1;
  printf("max: %d f: %d fshift %d ppm: %d Eb/No: %3.1f SNR: %3.1f bits: %d\r", max(s), f, shift, states.ppm, EbNo, SNR, length(rx_bits_buf));

  % look for complete Horus frame, delimited by 2 unique words

  bit = 1;
  nbits = length(rx_bits_buf);
  uw_loc1 = find_uw(states, bit, rx_bits_buf);

  if uw_loc1 != -1
    uw_loc2 = find_uw(states, uw_loc1+length(uw), rx_bits_buf);

    if uw_loc2 != -1

      % Now we can extract ascii chars from the frame

      str = [];
      st = uw_loc1 + length(states.uw);  % first bit of first char
      for i=st:nfield+npad:uw_loc2
        field = rx_bits_buf(i:i+nfield-1);
        ch_dec = field * (2.^(0:nfield-1))';
        % filter out unlikely characters that bit errors may introduce, and ignore \n
        if (ch_dec > 31) && (ch_dec < 91)
          str = [str char(ch_dec)];
        else 
          str = [str char(32)]; % space is "not sure"
        end
      end
      printf("%s\n", str);

      % throw out used bits in buffer

      rx_bits_buf =  rx_bits_buf(uw_loc2-1:length(rx_bits_buf));
    end
  end
  [s,c] = fread(stdin, N, "short");

endwhile
