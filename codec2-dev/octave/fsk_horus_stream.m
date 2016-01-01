#!/usr/bin/octave -qf

% fsk_horus_stream.m
% David Rowe 13 Oct 2015
%
% Experimental real time near space balloon FSK demodulator, takes
% 8kHz 16 bit samples from stdin, output txt string on stdout
%
% usage:
%  $ chmod 777 fsk_horus_stream.m
%  $ rec -t raw -r 8000 -s -2 -c 1 - -q | ./fsk_horus_stream.m
%  or (for those of us that avoid alsa like the plague)
%  $ arecord -D pulse -r 8000 -c 1 -f S16_LE - | ./fsk_horus_stream.m
%  and use the 'pavucontrol' utility to select a sound device for arecord.
%
%
% OR to test with a stored file (8kHz 16-bit shorts):
%  $ cat ~/Desktop/vk5arg-3.wav | ./fsk_horus_stream.m
%

% include modem library

fsk_horus_as_a_lib = 1; % make sure calls to test functions at bottom are disabled
fsk_horus;  

gps_log = "~/Desktop/gps_log.txt"
system_command = "echo -n \"/home/david/Desktop/gps_log.txt\" | nc -u -q1 127.0.0.1 21234";

% Upload Telemetry to Habitat (http://tracker.habhub.org/)
telem_upload_enabled = false;
% Update this command with your own callsign.
telem_upload_command = "python telem_upload.py -c N0CALL_Octave";

more off;
states = fsk_horus_init(8000, 100);
%states = fsk_horus_init_rtty_uw(states);
N = states.N;
Rs = states.Rs;
nsym = states.nsym;
nin = states.nin;
nfield = states.rtty.nfield;
npad = states.rtty.npad;
EbNo = 0;
SNR = 0;

rx = [];
rx_bits_buf = [];
rx_bits_sd_buf = [];

[s,c] = fread(stdin, N, "short");

while c

  rx = [rx s'];
 
  % demodulate samples to bit stream

  while length(rx) > nin
    [rx_bits states] = fsk_horus_demod(states, rx(1:nin)');
    rx_bits_buf = [rx_bits_buf rx_bits];
    rx_bits_sd_buf = [rx_bits_sd_buf states.rx_bits_sd];
    rx = rx(nin+1:length(rx));
    nin = states.nin;
    EbNo = 0.9*EbNo + 0.1*states.EbNodB;
    SNR = EbNo + 10*log10(states.Rs/3000);
    %printf("nin: %d length(rx): %d length(rx_bits_buf): %d \n", nin, length(rx), length(rx_bits_buf));
  endwhile
  f = (states.f1+states.f2)/2; shift = states.f2 - states.f1;
  printf("max: %d f: %d fshift %d ppm: %d Eb/No: %3.1f SNR: %3.1f bits: %d\r", max(s), f, shift, states.ppm, EbNo, SNR, length(rx_bits_buf));

  packet_found = 0;

  % Look for complete Horus RTTY frame -------------------------------------------------------------------

  nbits = length(rx_bits_buf);
  uw_loc = find_uw(states.rtty, 1, rx_bits_buf);

  if uw_loc != -1
    packet_found = 1;

    if (uw_loc + states.rtty.max_packet_len) < nbits
      %printf("\n%d nbits: %d\n",uw_loc + states.rtty.max_packet_len, nbits);

      [str crc_ok] = extract_ascii(states.rtty, rx_bits_buf, uw_loc);

      if crc_ok == 0
        [str_flipped crc_flipped_ok] = sd_bit_flipping(states.rtty, rx_bits_buf, rx_bits_sd_buf, uw_loc, uw_loc+states.rtty.max_packet_len); 
        if crc_flipped_ok
          str = sprintf("%s fixed", str_flipped);
          crc_ok = 1;
        end
      end
      
      if crc_ok
        if telem_upload_enabled
          % Upload to Habitat.
          ascii_upload_cmd = sprintf("%s %s",telem_upload_command,str);
          printf("Uploading ASCII to Habitat...\n");
          system(ascii_upload_cmd,false,"async");
        end

        strok = sprintf("%s CRC OK", str);
      else
        strok = sprintf("%s CRC BAD", str);
      end
      printf("\n  %s         \n", strok);
      
      % throw out used bits in buffer.  We're not sure where the next packet starts
      % so lets remove everything up to just after the UW we just used to force
      % a search for the next UW.

      rx_bits_buf    = rx_bits_buf(uw_loc+length(states.rtty.uw):length(rx_bits_buf));
      rx_bits_sd_buf = rx_bits_sd_buf(uw_loc+length(states.rtty.uw):length(rx_bits_sd_buf));

      if crc_ok
        % extract GPS coords and save to log file for mapping software

        str_split = strsplit(str,",");
        if length(str_split) > 4
          lat = str_split{1,4}; long = str_split{1,5};
          f = fopen(gps_log,"at");
          fprintf(f,"%s,%s\n", lat, long);
          fclose(f);
        end

        % TODO: thin out log file to max_points to lighten plotting load

        % tell foxtrotGPS to plot track

        system(system_command);
      end
    end
  end

  % Look for complete Horus BINARY frame -------------------------------------------------------------------

  nbits = length(rx_bits_buf);
  uw_loc = find_uw(states.binary, 1, rx_bits_buf);

  if uw_loc != -1
    packet_found = 1;

    if (uw_loc + states.binary.max_packet_len) < nbits

      pin = uw_loc;    
      for i=1:45
        rx_bytes(i) = rx_bits_buf(pin:pin+7) * (2.^(7:-1:0))';
        pin += 8;
        %printf("%d 0x%02x\n", i, rx_bytes(i));
      end

      printf("\n  ");
      f=fopen("horus_rx_bits_binary.txt","wt");
      fwrite(f, rx_bytes, "uchar");
      fclose(f);

      % horus_l2 can be compiled a bunch of different ways.  You need to
      % compile with:
      %   codec2-dev/src$ gcc horus_l2.c -o horus_l2 -Wall -DDEC_RX_BITS -DHORUS_L2_RX

      system("../src/horus_l2");
      if telem_upload_enabled
        % Upload binary payload data to Habitat.
        binary_upload_addition = "`cat horus_rx_bits_hex.txt`";
        binary_upload_cmd = sprintf("%s %s",telem_upload_command,binary_upload_addition);
        printf("Uploading Binary to Habitat...\n");
        system(binary_upload_cmd,type="async");
      end

      % throw out used bits in buffer.  We're not sure where the next packet starts
      % so lets remove everything up to just after the UW we just used to force
      % a search for the next UW.

      rx_bits_buf    = rx_bits_buf(uw_loc+length(states.binary.uw):length(rx_bits_buf));
      rx_bits_sd_buf = rx_bits_sd_buf(uw_loc+length(states.binary.uw):length(rx_bits_sd_buf));
    end
  end

  % Truncate buffers if no UW found so they don't grow endlessly with no signal.
  % Keep very end of it as it may have part of a UW in it

  if packet_found == 0
    max_len = states.rtty.max_packet_len*4;
    if length(rx_bits_buf) > max_len
      rx_bits_buf = rx_bits_buf(length(rx_bits_buf)-states.rtty.max_packet_len:length(rx_bits_buf));
      rx_bits_sd_buf = rx_bits_sd_buf(length(rx_bits_sd_buf)-states.rtty.max_packet_len:length(rx_bits_sd_buf));
    end
  end

  [s,c] = fread(stdin, N, "short");

endwhile
