% ofdm_tx.m
% David Rowe April 2017
%
% File based ofdm tx.  Generate a file of ofdm samples, inclduing
% optional channel simulation.


function ofdm_tx(filename, Nsec, EbNodB=100, channel='awgn', freq_offset_Hz=0)
  ofdm_lib;
  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 16; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  ofdm_load_const;

  % generate fixed test frame of tx bits and run OFDM modulator
  % todo: maybe extend this to 4 or 8 frames, one is a bit short

  Nrows = Nsec*Rs;
  Nframes = floor((Nrows-1)/Ns);

  rand('seed', 100);
  tx_bits = rand(1,Nbitsperframe) > 0.5;

  tx = [];
  for f=1:Nframes
    tx = [tx ofdm_mod(states, tx_bits)];
  end
  Nsam = length(tx);

  % channel simulation

  EbNo = bps * (10 .^ (EbNodB/10));
  variance = 1/(M*EbNo/2);
  woffset = 2*pi*freq_offset_Hz/Fs;

  % SNR calculation is probably a bit off or confused, need to compare
  % to cohpsk_frame_design.ods and make sure I've got it right.  PLus
  % 3 is used as we are dealing with real channels because mumble
  % mumble not understood very well mumble magic number.

  SNRdB = EbNodB + 10*log10(bps*Rs*Nc/Fs) + 3;
  printf("EbNo: %3.1f SNR(3k): %3.1f foff: %3.1f\n", EbNodB, SNRdB, freq_offset_Hz);

  % set up HF model ---------------------------------------------------------------

  if strcmp(channel, 'hf')

    % some typical values, or replace with user supplied

    dopplerSpreadHz = 1.0; path_delay_ms = 1;

    path_delay_samples = path_delay_ms*Fs/1000;
    printf("Doppler Spread: %3.2f Hz Path Delay: %3.2f ms %d samples\n", dopplerSpreadHz, path_delay_ms, path_delay_samples);

    % generate same fading pattern for every run

    randn('seed',1);

    spread1 = doppler_spread(dopplerSpreadHz, Fs, (Nsec*(M+Ncp)/M+0.2)*Fs);
    spread2 = doppler_spread(dopplerSpreadHz, Fs, (Nsec*(M+Ncp)/M+0.2)*Fs);

    % sometimes doppler_spread() doesn't return exactly the number of samples we need
 
    assert(length(spread1) >= Nsam, "not enough doppler spreading samples");
    assert(length(spread2) >= Nsam, "not enough doppler spreading samples");

    hf_gain = 1.0/sqrt(var(spread1)+var(spread2));
  end

  rx = tx;

  if strcmp(channel, 'hf')
    rx  = hf_gain * tx(1:Nsam) .* spread1(1:Nsam);
    rx += hf_gain * [zeros(1,path_delay_samples) tx(1:Nsam-path_delay_samples)] .* spread2(1:Nsam);
  end

  rx = rx .* exp(j*woffset*(1:Nsam));

  % note variance/2 as we are using real() operator

  noise = sqrt(variance/2)*(0.5*randn(1,Nsam) + j*0.5*randn(1,Nsam));
  rx += noise;
  10*log10(var(tx)/var(noise))

  Ascale = 4E5;
  frx=fopen(filename,"wb"); fwrite(frx, Ascale*real(rx), "short"); fclose(frx);
endfunction
