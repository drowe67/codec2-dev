% fdmdv_mod.m
%
% Modulator function for FDMDV modem, uses test frames as input and
% outputs a raw file of 16 bit shorts at a sample rate of 8 kHz.
%
% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

function tx_fdm = fdmdv_mod(rawfilename, nbits)

  fdmdv; % include modem code
  f = fdmdv_init;
  Nc = f.Nc; Nb = f.Nb;
  
  frames = floor(nbits/(Nc*Nb))
  tx_fdm = [];
  gain = 1000; % Scale up to 16 bit shorts
  prev_tx_symbols = ones(Nc+1,1); prev_tx_symbols(Nc+1) = 2;

  for i=1:frames
    [tx_bits f] = get_test_bits(f,Nc*Nb);
    [tx_symbols f] = bits_to_psk(f,prev_tx_symbols, tx_bits);
    prev_tx_symbols = tx_symbols;
    [tx_baseband f] = tx_filter(f, tx_symbols);
    tx_fdm = [tx_fdm real(fdm_upconvert(f, tx_baseband))];
  end

  tx_fdm *= gain;
  fout = fopen(rawfilename,"wb");
  fwrite(fout, tx_fdm, "short");
  fclose(fout);
endfunction
