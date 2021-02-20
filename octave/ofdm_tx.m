% ofdm_tx.m
% David Rowe March 2018
%
% File based, uncoded OFDM tx.  Generates a file of ofdm samples,
% including optional channel simulation.  See also ofdm_ldpc_tx.m, and
% ofdm_mod.c

#{
  Examples:

  i) 10 seconds, AWGN channel at SNR3k=3dB

    octave:4> ofdm_tx("awgn_snr_3dB_700d.raw", "700D", 10, 3);

  ii) 10 seconds, multipath poor channel at SNR=6dB

    octave:5> ofdm_tx("hf_snr_6dB_700d.raw", "700D", 10, 6, "mpp");
#}

function ofdm_tx(filename, mode="700D", Nsec, SNR3kdB=100, channel='awgn', freq_offset_Hz=0, tx_clip_en=0)
  ofdm_lib;
  channel_lib;
  randn('seed',1);
  pkg load signal;

  dpsk = 0;
  if strcmp(mode,"700D-DPSK")
    mode = "700D"; dpsk = 1;
  end
  if strcmp(mode,"2020-DPSK")
    mode = "2020"; dpsk = 1;
  end

  % init modem

  config = ofdm_init_mode(mode);
  states = ofdm_init(config);
  print_config(states);
  ofdm_load_const;
  states.dpsk=dpsk;

  % Generate fixed test frame of tx bits and run OFDM modulator

  Npackets = round(Nsec/states.Tpacket);
  tx_bits = create_ldpc_test_frame(states, coded_frame=0);
  atx = ofdm_mod(states, tx_bits);
  tx = [];
  for f=1:Npackets
    tx = [tx atx];
  end
  if states.data_mode
    tx = [states.tx_preamble tx];
  end

  printf("Npackets: %d  ", Npackets);
  states.verbose=1;
  tx = ofdm_hilbert_clipper(states, tx, tx_clip_en);
  rx_real = ofdm_channel(states, tx, SNR3kdB, channel, freq_offset_Hz);
  frx = fopen(filename,"wb"); fwrite(frx, rx_real, "short"); fclose(frx);
endfunction
