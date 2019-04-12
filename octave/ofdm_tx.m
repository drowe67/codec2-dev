% ofdm_tx.m
% David Rowe March 2018
%
% File based, uncoded OFDM tx.  Generates a file of ofdm samples,
% including optional channel simulation.  See also ofdm_ldpc_tx.m, and
% ofdm_mod.c

#{
  Examples:
 
  i) 10 seconds, AWGN channel at Eb/No=3dB

    octave:4> ofdm_tx('awgn_ebno_3dB_700d.raw', "700D", 10, 3);

  ii) 10 seconds, HF channel at Eb/No=6dB

    ofdm_tx('hf_ebno_6dB_700d.raw',  "700D", 10, 6, 'hf');
    
  iii) 10 seconds, 2200 waveform, AWGN channel, Eb/No=100dB (effectively noise free)

    ofdm_tx('hf_ebno_6dB_700d.raw',  "2200", 10);
#}

% Note EbNodB is for payload data bits, so will be 10log10(rate) higher than
% raw EbNodB used in ofdm_tx() at uncoded bit rate

function ofdm_tx(filename, mode="700D", Nsec, EbNodB=100, channel='awgn', freq_offset_Hz=0, dfoff_hz_per_sec = 0, initial_noise_sams=0, tx_filter=0)
  ofdm_lib;

  % init modem
  
  [bps Rs Tcp Ns Nc] = ofdm_init_mode(mode)
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  print_config(states);
  ofdm_load_const;

  % Generate fixed test frame of tx bits and run OFDM modulator

  Nrows = Nsec*Rs;
  Nframes = floor((Nrows-1)/Ns);
  tx_bits = create_ldpc_test_frame(states, coded_frame=0);

  tx = [];
  for f=1:Nframes
    tx = [tx ofdm_mod(states, tx_bits)];
  end
  
  Nsam = length(tx);

  % channel simulation

  EsNo = rate * bps * (10 .^ (EbNodB/10));
  variance = 1/(M*EsNo/2);
  woffset = 2*pi*freq_offset_Hz/Fs;
  dwoffset = 2*pi*dfoff_hz_per_sec/(Fs*Fs);
  
  SNRdB = EbNodB + 10*log10(Nc*bps*Rs/3000);
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

  % experimental coarse amplitude quantisation

  quant_tx = 0;
  if quant_tx
    tx_re = real(tx); tx_im = imag(tx);
    tx_re = min(tx_re,0.5); tx_re = max(tx_re,-0.5);
    tx_im = min(tx_im,0.5); tx_im = max(tx_im,-0.5);
    step = 0.05/4;
    tx_re = step*round(tx_re/step);
    tx_im = step*round(tx_im/step);
    tx = tx_re + j*tx_im;
    figure(1); clf; subplot(211); plot(real(tx(1:100)),'+-'); subplot(212); plot(imag(tx(1:100)),'+-'); 
  end
  
  rx = tx;
  if tx_filter
    bpf_coeff = make_ofdm_bpf(write_c_header_file=0);
    rx = filter(bpf_coeff,1,tx);
    figure(1); clf;
    subplot(211);
    plot((1:length(tx))*8000/length(tx), 20*log10(abs(fft(tx))))
    axis([1 4000 0 60])
    subplot(212);
    plot((1:length(tx))*8000/length(tx), 20*log10(abs(fft(rx))))
    axis([1 4000 0 60])
  end
  
  if strcmp(channel, 'hf')
    rx  = tx(1:Nsam) .* spread1(1:Nsam);
    rx += [zeros(1,path_delay_samples) tx(1:Nsam-path_delay_samples)] .* spread2(1:Nsam);

    % normalise rx power to same as tx

    nom_rx_pwr = 2/(Ns*(M*M)) + Nc/(M*M);
    rx_pwr = var(rx);
    rx *= sqrt(nom_rx_pwr/rx_pwr);
  end

  phase_offset = woffset*(1:Nsam) + 0.5*dwoffset*((1:Nsam).^2);
  rx = rx .* exp(j*phase_offset);

  rx = [zeros(1,initial_noise_sams) rx];
  Nsam = length(rx);
  
  % note variance/2 as we are using real() operator, mumble,
  % reflection of -ve freq to +ve, mumble, hand wave

  randn('seed',1);
  noise = sqrt(variance/2)*0.5*randn(1,Nsam);
  rx = real(rx) + noise;
  printf("measured SNR: %3.2f dB\n", 10*log10(var(real(tx))/var(noise))+10*log10(4000) - 10*log10(3000));

  % adjusted by experiment to match rms power of early test signals

  frx=fopen(filename,"wb"); fwrite(frx, states.amp_scale*rx, "short"); fclose(frx);
endfunction
