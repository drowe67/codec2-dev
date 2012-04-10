% fdmdv_demod.m
%
% Demodulator function for FDMDV modem.  Requires 48kHz sample rate raw files
% as input
%
% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

function fdmdv_demod(rawfilename, nbits)

  fdmdv; % include modem code

  modulation = 'dqpsk';

  fin = fopen(rawfilename, "rb");
  gain = 1000;
  frames = nbits/(Nc*Nb);

  prev_rx_symbols = ones(Nc+1,1);
  foff_phase = 1;

  % pilot LUT, used for copy of pilot at rx
  
  pilot_lut = generate_pilot_lut;
  pilot_lut_index = 1;
  prev_pilot = zeros(1,M);

  % BER stats

  total_bit_errors = 0;
  total_bits = 0;
  bit_errors_log = [];
  sync_log = [];
  test_frame_sync_log = [];
  test_frame_sync_state = 0;

  rx_symbols_log = [];
  rx_timing_log = [];
  foff_log = [];

  % resampler states

  t = 3;
  ratio = 1.000;
  F=6;
  MF=M*F;
  nin = MF;
  nin_size = MF+6;
  buf_in = zeros(1,nin_size);
  rx_fdm_buf = [];
  ratio_log = [];

  % Main loop ----------------------------------------------------

  for f=1:frames
    % update buf_in memory

    m = nin_size - nin;
    for i=1:m
      buf_in(i) = buf_in(i+nin);  
    end
    
    % obtain nin samples of the test input signal

    for i=m+1:nin_size
      buf_in(i) = fread(fin, 1, "short")/gain; 
    end

    % resample at 48kHz and decimate to 8kHz

    [rx_fdm_mf t nin] = resample(buf_in, t, 1.0, MF);
    rx_fdm = rx_fdm_mf(1:F:MF);

    %rx_fdm = buf_in(m+1:m+n);

    %for i=1:M
    %  rx_fdm(i) = fread(fin, 1, "short")/gain; 
    %end
    rx_fdm_buf = [rx_fdm_buf rx_fdm];

    % frequency offset estimation and correction

    for i=1:M
      pilot(i) = pilot_lut(pilot_lut_index);
      pilot_lut_index++;
      if pilot_lut_index > 4*M
        pilot_lut_index = 1;
      end
    end
    foff = rx_est_freq_offset(rx_fdm, pilot, prev_pilot);
    prev_pilot = pilot;
    foff_log = [ foff_log foff ];
    foff = 0;
    foff_rect = exp(j*2*pi*foff/Fs);

    for i=1:M
      foff_phase *= foff_rect';
      rx_fdm(i) = rx_fdm(i)*foff_phase;
    end

    % baseband processing

    rx_baseband = fdm_downconvert(rx_fdm);
    rx_filt = rx_filter(rx_baseband);

    [rx_symbols rx_timing] = rx_est_timing(rx_filt, rx_baseband);
    rx_timing_log = [rx_timing_log rx_timing];
    beta = 1E-7;
    ratio += beta*rx_timing;
    if (ratio > 1.002)
       ratio = 1.002;
    end
    if (ratio < 0.998)
       ratio = 0.998;
    end
    ratio_log = [ratio_log ratio];

    if strcmp(modulation,'dqpsk')
      rx_symbols_log = [rx_symbols_log rx_symbols.*conj(prev_rx_symbols)*exp(j*pi/4)];
    else
      rx_symbols_log = [rx_symbols_log rx_symbols];
    endif
    [rx_bits sync] = qpsk_to_bits(prev_rx_symbols, rx_symbols, modulation);
    prev_rx_symbols = rx_symbols;
    sync_log = [sync_log sync];

    % count bit errors if we find a test frame
    % Allow 15 frames for filter memories to fill and time est to settle

    [test_frame_sync bit_errors] = put_test_bits(rx_bits);
    if (test_frame_sync == 1 && f > 15)
      total_bit_errors = total_bit_errors + bit_errors;
      total_bits = total_bits + Ntest_bits;
      bit_errors_log = [bit_errors_log bit_errors];
    end

    % test frame sync state machine, just for more informative plots
    
    next_test_frame_sync_state = test_frame_sync_state;
    if (test_frame_sync_state == 0)
      if (test_frame_sync == 1)      
        next_test_frame_sync_state = 1;
	test_frame_count = 0;
      end
    end

    if (test_frame_sync_state == 1)
      % we only expect another test_frame_sync pulse every 4 symbols
      test_frame_count++;
      if (test_frame_count == 4)
        test_frame_count = 0;
        if ((test_frame_sync == 0))      
          next_test_frame_sync_state = 0;
        end
      end
    end
    test_frame_sync_state = next_test_frame_sync_state;
    test_frame_sync_log = [test_frame_sync_log test_frame_sync_state];

  end

  % ---------------------------------------------------------------------
  % Print Stats
  % ---------------------------------------------------------------------

  ber = total_bit_errors / total_bits;

  printf("%d bits  %d errors  BER: %1.4f\n",total_bits, total_bit_errors, ber);

  % ---------------------------------------------------------------------
  % Plots
  % ---------------------------------------------------------------------

  figure(1)
  clf;
  [n m] = size(rx_symbols_log);
  plot(real(rx_symbols_log(1:Nc+1,10:m)),imag(rx_symbols_log(1:Nc+1,10:m)),'+')
  axis([-2 2 -2 2]);
  title('Scatter Diagram');

  figure(2)
  clf;
  subplot(211)
  plot(rx_timing_log(15:m))
  title('timing offset (samples)');
  subplot(212)
  plot(foff_log)
  title('Freq offset (Hz)');
  grid

  figure(3)
  clf;
  subplot(211)
  %plot(rx_fdm_buf);
  %title('FDM Rx Signal');
  plot(ratio_log-1);
  title('Sampling Clock error (ppm)');
  subplot(212)
  Nfft=Fs;
  S=fft(rx_fdm_buf,Nfft);
  SdB=20*log10(abs(S));
  plot(SdB(1:Fs/4))
  title('FDM Tx Spectrum');

  figure(4)
  clf;
  subplot(311)
  stem(sync_log)
  axis([0 frames 0 1.5]);
  title('BPSK Sync')
  subplot(312)
  stem(bit_errors_log);
  title('Bit Errors for test frames')
  subplot(313)
  plot(test_frame_sync_log);
  axis([0 frames 0 1.5]);
  title('Test Frame Sync')

  mean(ratio_log)
endfunction
