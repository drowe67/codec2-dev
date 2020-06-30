% ofdm_ldpc_tx.m
% David Rowe April 2017
%
% File based ofdm tx with LDPC encoding and interleaver.  Generates a
% file of ofdm samples, including optional channel simulation.

function ofdm_ldpc_tx(filename, mode="700D", Nsec, SNR3kdB=100, channel='awgn', freq_offset_Hz=0)
  ofdm_lib;
  ldpc;
  gp_interleaver;
  channel_lib;
  randn('seed',1);
  more off;
  
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
  elseif strcmp(mode, "datac1") || strcmp(mode, "datac2") || strcmp(mode, "datac3") || strcmp(mode, "qam16")
    payload_bits = round(ofdm_rand(code_param.data_bits_per_frame)/32767);
  end
  [packet_bits bits_per_packet] = fec_encode(states, code_param, mode, payload_bits, Ncodecframespermodemframe, Nbitspercodecframe);
   
  % modulate to create symbols and interleave  
  tx_bits = tx_symbols = [];
  tx_bits = [tx_bits payload_bits];
  for b=1:2:bits_per_packet
    tx_symbols = [tx_symbols qpsk_mod(packet_bits(b:b+1))];
  end
  assert(gp_deinterleave(gp_interleave(tx_symbols)) == tx_symbols);
  tx_symbols = gp_interleave(tx_symbols);
  
  % generate txt (non FEC protected) symbols
  txt_bits = zeros(1,Ntxtbits);
  txt_symbols = [];
  for b=1:2:length(txt_bits)
    txt_symbols = [txt_symbols qpsk_mod(txt_bits(b:b+1))];
  end

  % assemble interleaved modem packet that include UW and txt symbols  
  modem_packet = assemble_modem_packet_symbols(states, tx_symbols, txt_symbols);
  atx = ofdm_txframe(states, modem_packet); tx = [];
  for f=1:Npackets
    tx = [tx atx];
  end
  Nsam = length(tx);

  printf("Packets: %3d SNR(3k): %3.1f dB foff: %3.1f Hz ", Npackets, SNR3kdB, freq_offset_Hz);
  rx = channel_simulate(Fs, SNR3kdB, freq_offset_Hz, channel, tx);
  frx=fopen(filename,"wb"); fwrite(frx, states.amp_scale*rx, "short"); fclose(frx);
endfunction
