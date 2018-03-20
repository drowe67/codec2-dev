% ofdm_tx.m
% David Rowe March 2018
%
% File based ofdm tx.  Generates a file of ofdm samples, including
% optional channel simulation.  See also ofdm_ldpc_tx.m

#{
  Examples:
 
  i) 10 seconds, AWGN channel at Eb/No=3dB

    octave:4> ofdm_tx('awgn_ebno_3dB_700d.raw', 10,3);

  ii) 10 seconds, HF channel at Eb/No=6dB

    ofdm_tx('hf_ebno_6dB_700d.raw', 10, 6, 'hf');
#}


function ofdm_tx(filename, Nsec, EbNodB=100, channel='awgn', freq_offset_Hz=0)
  ofdm_lib;

  % init modem

  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 16; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;

  % Generate fixed test frame of tx bits and run OFDM modulator

  Nrows = Nsec*Rs;
  Nframes = floor((Nrows-1)/Ns);
  rand('seed', 100);
  tx_bits = round(rand(1,Nbitsperframe));

  tx = [];
  for f=1:Nframes
    tx = [tx ofdm_mod(states, tx_bits)];
  end
  
  Nsam = length(tx);

  % channel simulation

  EsNo = rate * bps * (10 .^ (EbNodB/10));
  variance = 1/(M*EsNo/2);
  woffset = 2*pi*freq_offset_Hz/Fs;

  SNRdB = EbNodB + 10*log10(700/3000);
  printf("EbNo: %3.1f dB  SNR(3k) est: %3.1f dB  foff: %3.1fHz ", EbNodB, SNRdB, freq_offset_Hz);

  % set up HF model ---------------------------------------------------------------

  if strcmp(channel, 'hf')
    randn('seed',1);

    % some typical values, or replace with user supplied

    dopplerSpreadHz = 1; path_delay_ms = 1;

    path_delay_samples = path_delay_ms*Fs/1000;
    printf("Doppler Spread: %3.2f Hz Path Delay: %3.2f ms %d samples\n",
           dopplerSpreadHz, path_delay_ms, path_delay_samples);

    % generate same fading pattern for every run

    randn('seed',1);

    spread1 = doppler_spread(dopplerSpreadHz, Fs, (Nsec*(M+Ncp)/M)*Fs*1.1);
    spread2 = doppler_spread(dopplerSpreadHz, Fs, (Nsec*(M+Ncp)/M)*Fs*1.1);

    % sometimes doppler_spread() doesn't return exactly the number of samples we need
 
    assert(length(spread1) >= Nsam, "not enough doppler spreading samples");
    assert(length(spread2) >= Nsam, "not enough doppler spreading samples");
  end

  rx = tx;

  if strcmp(channel, 'hf')
    rx  = tx(1:Nsam) .* spread1(1:Nsam);
    rx += [zeros(1,path_delay_samples) tx(1:Nsam-path_delay_samples)] .* spread2(1:Nsam);

    % normalise rx power to same as tx

    nom_rx_pwr = 2/(Ns*(M*M)) + Nc/(M*M);
    rx_pwr = var(rx);
    rx *= sqrt(nom_rx_pwr/rx_pwr);
  end

  rx = rx .* exp(j*woffset*(1:Nsam));

  % note variance/2 as we are using real() operator, mumble,
  % reflection of -ve freq to +ve, mumble, hand wave

  noise = sqrt(variance/2)*0.5*randn(1,Nsam);
  rx = real(rx) + noise;
  printf("measured SNR: %3.2f dB\n", 10*log10(var(real(tx))/var(noise))+10*log10(4000) - 10*log10(3000));

  % adjusted by experiment to match rms power of early test signals

  Ascale = 2E5*1.1491;

  frx=fopen(filename,"wb"); fwrite(frx, Ascale*rx, "short"); fclose(frx);
endfunction
