% ofdm_tx.m
% David Rowe April 2017
%
% File based ofdm tx.  Generate a file of ofdm samples, inclduing
% optional channel simulation.


function ofdm_tx(filename, Nsec)
  ofdm_lib;
  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 16; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;

  % generate fixed test frame of tx bits and run OFDM modulator

  Nrows = Nsec*Rs;
  Nframes = floor((Nrows-1)/Ns);

  rand('seed', 100);
  tx_bits = rand(1,Nbitsperframe) > 0.5;

  tx = [];
  for f=1:Nframes
    tx = [tx ofdm_mod(states, tx_bits)];
  end

  Ascale = 32000/max(real(tx));
  ftx=fopen(filename,"wb"); fwrite(ftx, Ascale*real(tx), "short"); fclose(ftx);
endfunction
