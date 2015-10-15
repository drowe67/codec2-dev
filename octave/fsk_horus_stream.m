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

% include library (make sure calls to test functions at bottom are disabled)

fsk_horus_as_a_lib = 1;
fsk_horus;  

gps_log = "~/Desktop/gps_log.txt"
system_command = "echo -n \"/home/david/Desktop/gps_log.txt\" | nc -u -q1 127.0.0.1 21234"

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
      [str crc_ok] = extract_ascii(states, rx_bits_buf, uw_loc1, uw_loc2);
      printf("%s         \n", str);

      % throw out used bits in buffer

      rx_bits_buf =  rx_bits_buf(uw_loc2-1:length(rx_bits_buf));

      % extract GPS coords and save to log file for mapping software

      str_split = strsplit(str,",");
      lat = str_split{1,4}; long = str_split{1,5};
      f = fopen(gps_log,"at");
      fprintf(f,"%s,%s\n", lat, long);
      fclose(f);

      % thin out log file to points_max points to lighten plotting load
      % TBC

      % tell foxtrotGPS to plot track

      system(system_command);
    end
  end
  [s,c] = fread(stdin, N, "short");

endwhile
