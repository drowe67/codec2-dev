% ofdm_ldpc_tx.m
% David Rowe April 2017
%
% File based ofdm tx with LDPC encoding and interleaver.  Generates a
% file of ofdm samples, including optional channel simulation.

function ofdm_ldpc_tx(filename, mode="700D", Nsec, SNR3kdB=100, channel='awgn', freq_offset_Hz=0)
  ofdm_lib;
  ldpc;
  gp_interleaver;

  % init modem

  config = ofdm_init_mode(mode);
  states = ofdm_init(config);
  print_config(states);
  ofdm_load_const;

  % some constants used for assembling modem frames
  
  [code_param Nbitspercodecframe Ncodecframespermodemframe] = codec_to_frame_packing(states, mode);

  % Generate fixed test frame of tx bits and run OFDM modulator

  Npackets = round(Nsec/states.Tpacket);

  % OK generate a modem frame using random payload bits

  if strcmp(mode, "700D")
    payload_bits = round(ofdm_rand(code_param.data_bits_per_frame)/32767);
  elseif strcmp(mode, "2020")
    payload_bits = round(ofdm_rand(Ncodecframespermodemframe*Nbitspercodecframe)/32767);
  elseif strcmp(mode, "datac1") || strcmp(mode, "datac2") || strcmp(mode, "datac3")
    payload_bits = round(ofdm_rand(code_param.data_bits_per_frame)/32767);
  end
  [frame_bits bits_per_frame] = assemble_frame(states, code_param, mode, payload_bits, Ncodecframespermodemframe, Nbitspercodecframe);
   
  % modulate to create symbols and interleave
  
  tx_bits = tx_symbols = [];
  tx_bits = [tx_bits payload_bits];
  for b=1:2:bits_per_frame
    tx_symbols = [tx_symbols qpsk_mod(frame_bits(b:b+1))];
  end
  assert(gp_deinterleave(gp_interleave(tx_symbols)) == tx_symbols);
  tx_symbols = gp_interleave(tx_symbols);
  
  % generate txt symbols
 
  txt_bits = zeros(1,Ntxtbits);
  txt_symbols = [];
  for b=1:2:length(txt_bits)
    txt_symbols = [txt_symbols qpsk_mod(txt_bits(b:b+1))];
  end

  % assemble interleaved modem packet that include UW and txt symbols
  
  modem_frame = assemble_modem_packet_symbols(states, tx_symbols, txt_symbols);
  atx = ofdm_txframe(states, modem_frame); tx = [];
  for f=1:Npackets
    tx = [tx atx];
  end
  % a few empty frames of samples os Rx can finish it's processing
  tx = [tx zeros(1,2*Nsamperframe)]; 
  Nsam = length(tx);

  printf("Packets: %3d SNR(3k): %3.1f dB foff: %3.1f Hz ", Npackets, SNR3kdB, freq_offset_Hz);

  % set up HF model ---------------------------------------------------------------

  if strcmp(channel, 'awgn') == 0
    randn('seed',1);

    % Winlink multipath definitions
    if strcmp(channel, 'mpg')     dopplerSpreadHz = 0.1; path_delay_ms = 0.5;
    elseif strcmp(channel, 'mpm') dopplerSpreadHz = 0.5; path_delay_ms = 1.0;
    elseif strcmp(channel, 'mpp') dopplerSpreadHz = 1.0; path_delay_ms = 2.0;
    elseif strcmp(channel, 'mpd') dopplerSpreadHz = 2.5; path_delay_ms = 4.0;
    elseif printf("Unknown multipath channel\n"); assert(0); end
    
    path_delay_samples = path_delay_ms*Fs/1000;
    printf(" Doppler Spread: %3.2f Hz Path Delay: %3.2f ms %d samples\n", dopplerSpreadHz, path_delay_ms, path_delay_samples);

    % generate same fading pattern for every run
    randn('seed',1);

    spread1 = doppler_spread(dopplerSpreadHz, Fs, Fs*Nsec + 2*Np*Nsamperframe);
    spread2 = doppler_spread(dopplerSpreadHz, Fs, Fs*Nsec + 2*Np*Nsamperframe);
   
    % sometimes doppler_spread() doesn't return exactly the number of samples we need 
    if length(spread1) < Nsam
      printf("not enough doppler spreading samples %d %d\n", length(spread1), Nsam);
      assert(0);
    end
    if length(spread2) < Nsam
      printf("not enough doppler spreading samples %d %d\n", length(spread2), Nsam);
      assert(0);
    end
  end

  rx = tx;

  if strcmp(channel, 'awgn') == 0
    % multipath model, this will affect signal power but we take care of that in the SNR calculations below
    rx  = tx(1:Nsam) .* spread1(1:Nsam);
    rx += [zeros(1,path_delay_samples) tx(1:Nsam-path_delay_samples)] .* spread2(1:Nsam);
  end

  woffset = 2*pi*freq_offset_Hz/Fs;
  rx = rx .* exp(j*woffset*(1:Nsam));

  rx = real(rx); S = rx*rx';

  % SNR in a 4k bandwidth will be lower than 3k as total noise power N is higher
  SNR4kdB = SNR3kdB - 10*log10(Fs/2) + 10*log10(3000); SNR = 10^(SNR4kdB/10);
  N = S/SNR; sigma = sqrt(N/Nsam);
  n = sigma*randn(1,Nsam);
  % printf("SNR3kdB: %f SNR4kdB: %f N: %f %f\n", SNR3kdB, SNR4kdB, N, n*n');
  rx += n;
  % check our sums are OK to within 0.25 dB
  SNR4kdB_measured = 10*log10(S/(n*n')); assert (abs(SNR4kdB - SNR4kdB_measured) < 0.25);
  printf("meas SNR3k: %3.2f dB\n", 10*log10(S/(n*n')) + 10*log10(4000) - 10*log10(3000));

  frx=fopen(filename,"wb"); fwrite(frx, states.amp_scale*rx, "short"); fclose(frx);
endfunction
