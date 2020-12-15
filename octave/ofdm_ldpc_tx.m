% ofdm_ldpc_tx.m
% David Rowe April 2017
%
% File based ofdm tx with LDPC encoding and interleaver.  Generates a
% file of ofdm samples, including optional channel simulation.

function ofdm_ldpc_tx(filename, mode="700D", Nsec, SNR3kdB=100, channel='awgn', freq_offset_Hz=0, tx_clip_en=0)
  ofdm_lib;
  ldpc;
  gp_interleaver;
  channel_lib;
  pkg load signal;
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

  if strcmp(mode, "2020")
    payload_bits = round(ofdm_rand(Ncodecframespermodemframe*Nbitspercodecframe)/32767);
  else
    payload_bits = round(ofdm_rand(code_param.data_bits_per_frame)/32767);
  end
  [packet_bits bits_per_packet] = fec_encode(states, code_param, mode, payload_bits, Ncodecframespermodemframe, Nbitspercodecframe);

  % modulate to create symbols and interleave
  tx_symbols = [];
  for b=1:bps:bits_per_packet
    if bps == 2 tx_symbols = [tx_symbols qpsk_mod(packet_bits(b:b+bps-1))]; end
    if bps == 4 tx_symbols = [tx_symbols qam16_mod(states.qam16, packet_bits(b:b+bps-1))]; end
  end
  assert(gp_deinterleave(gp_interleave(tx_symbols)) == tx_symbols);
  tx_symbols = gp_interleave(tx_symbols);

  % generate txt (non FEC protected) symbols
  txt_bits = zeros(1,Ntxtbits);
  txt_symbols = [];
  for b=1:bps:length(txt_bits)
    if bps == 2 txt_symbols = [txt_symbols qpsk_mod(txt_bits(b:b+bps-1))]; end
    if bps == 4 txt_symbols = [txt_symbols qam16_mod(states.qam16,txt_bits(b:b+bps-1))]; end
  end

  % assemble interleaved modem packet that include UW and txt symbols
  modem_packet = assemble_modem_packet_symbols(states, tx_symbols, txt_symbols);

  % sanity check
  [rx_uw rx_codeword_syms payload_amps txt_bits] = disassemble_modem_packet(states, modem_packet, ones(1,length(modem_packet)));
  assert(rx_uw == states.tx_uw);

  atx = ofdm_txframe(states, modem_packet); tx = [];
  for f=1:Npackets
    tx = [tx atx];
  end

  printf("Npackets: %d  ", Npackets);
  [rx_real rx] = ofdm_clip_channel(states, tx, SNR3kdB, channel, freq_offset_Hz, tx_clip_en);
  frx=fopen(filename,"wb"); fwrite(frx, rx_real, "short"); fclose(frx);
  if length(rx) >= states.Fs
    figure(1); clf; plot(20*log10(abs(fft(rx(1:states.Fs)/16384))));
    axis([1 states.Fs -20 60])
  end
endfunction
