% fsk_demod_file.m
% David Rowe May 2020
%
% Demodulate a file of off air samples and plot a bunch of internal
% states. Useful for debugging the FSK demod configuration

function fsk_demod_file(filename, format="s16", Fs=8000, Rs=50, M=2, max_frames=1E32)
  more off;
  fsk_lib;
  fin = fopen(filename,"rb"); 
  read_complex = 0; sample_size = 'int16'; plot_en = 1;
  if strcmp(format,"cs16") read_complex = 1; end

  states = fsk_init(Fs, Rs, M);
  nbit = states.nbit;

  frames = 0;
  rx = []; rx_bits_log = []; rx_bits_sd_log = []; norm_rx_timing_log = [];
  f_int_resample_log = []; EbNodB_log = []; ppm_log = [];
  f_log = []; Sf_log = [];
  
  % Extract raw bits from samples ------------------------------------------------------

  printf("demod of raw bits....\n");

  finished = 0; ph = 1;
  while (finished == 0)

    % read nin samples from input file

    nin = states.nin;
    if read_complex
      [sf count] = fread(fin, 2*nin, sample_size);
      if sample_size == "uint8" sf = (sf - 127)/128; end
      sf = sf(1:2:end) + j*sf(2:2:end);
      count /= 2;
      if shift_fs_on_4
        % optional shift up in freq by Fs/4 to get into freq est range
        for i=1:count
          ph = ph*exp(j*pi/4);
          sf(i) *= ph;
        end
      end
    else
      [sf count] = fread(fin, nin, "short");
    end
    rx = [rx; sf];
    
    if count == nin
      frames++;

      % demodulate to stream of bits

      states = est_freq(states, sf, states.M);
      if states.freq_est_type == 'mask' states.f = states.f2; end
      [rx_bits states] = fsk_demod(states, sf);

      rx_bits_log = [rx_bits_log rx_bits];
      rx_bits_sd_log = [rx_bits_sd_log states.rx_bits_sd];
      norm_rx_timing_log = [norm_rx_timing_log states.norm_rx_timing];
      f_int_resample_log = [f_int_resample_log abs(states.f_int_resample)];
      EbNodB_log = [EbNodB_log states.EbNodB];
      ppm_log = [ppm_log states.ppm];
      f_log = [f_log; states.f];
      Sf_log = [Sf_log; states.Sf'];
    else
      finished = 1;
    end

    if frames > max_frames finished=1; end
      
  end
  printf("frames: %d\n", frames);
  fclose(fin);

  if plot_en
    printf("plotting...\n");

    figure(1); clf;
    plot(f_log);
    title('Tone Freq Estimates');
    
    figure(2);
    plot(f_int_resample_log','+')
    title('Integrator outputs for each tone');

    figure(3); clf
    subplot(211)
    plot(norm_rx_timing_log)
    axis([1 frames -0.5 0.5])
    title('norm fine timing')
    subplot(212)
    plot(states.nerr_log)
    title('num bit errors each frame')
 
    figure(4); clf
    plot(EbNodB_log);
    title('Eb/No estimate')

    figure(5); clf
    rx_nowave = rx(1000:length(rx)); % skip past wav header if it's a wave file
    subplot(211)
    plot(real(rx_nowave));
    title('input signal to demod (1 sec)')
    xlabel('Time (samples)');
    %axis([1 states.Fs -35000 35000])

    % normalise spectrum to 0dB full scale with sine wave input
    subplot(212);
    if sample_size == "int16" max_value = 32767; end
    if sample_size == "uint8" max_value = 127; end
    RxdBFS = 20*log10(abs(fft(rx_nowave(1:states.Fs)))) - 20*log10((states.Fs/2)*max_value);
    plot(RxdBFS)
    axis([1 states.Fs/2 -80 0])
    xlabel('Frequency (Hz)');

    figure(6); clf
    plot(ppm_log)
    title('Sample clock (baud rate) offset in PPM');

    figure(7); clf; mesh(Sf_log(1:10,:));
  end

endfunction

