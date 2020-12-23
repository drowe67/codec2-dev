% ofdm_tx.m
% David Rowe March 2018
%
% File based, uncoded OFDM tx.  Generates a file of ofdm samples,
% including optional channel simulation.  See also ofdm_ldpc_tx.m, and
% ofdm_mod.c

#{
  Examples:

  i) 10 seconds, AWGN channel at SNR3k=3dB

    octave:4> ofdm_tx("awgn_ebno_3dB_700d.raw", "700D", 10, 3);

  ii) 10 seconds, multipath poor channel at SNR=6dB

    ofdm_tx("hf_ebno_6dB_700d.raw", "700D", 10, 6, "mpp");

  iii) 10 seconds, 2200 waveform, AWGN channel, SNR3k=100dB (noise free)

    ofdm_tx("hf_2020.raw", "2200", 10);
#}

% Note EbNodB is for payload data bits, so will be 10log10(rate) higher than
% raw EbNodB used in ofdm_tx() at uncoded bit rate

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
  tx = [];
  for f=1:Npackets
    tx = [tx ofdm_mod(states, tx_bits)];
  end

  printf("Npackets: %d  ", Npackets);
  rx = ofdm_clip_channel(states, tx, SNR3kdB, channel, freq_offset_Hz, tx_clip_en);
  frx=fopen(filename,"wb"); fwrite(frx, rx, "short"); fclose(frx);
endfunction
