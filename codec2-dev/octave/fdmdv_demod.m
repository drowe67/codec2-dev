% fdmdv_demod.m
%
% Demodulator function for FDMDV modem (Octave version).  Requires
% 8kHz sample rate raw files as input
%
% Copyright David Rowe 2012
% This program is distributed under the terms of the GNU General Public License 
% Version 2
%

function fdmdv_demod(rawfilename, nbits, pngname)

  fdmdv; % include modem code

  modulation = 'dqpsk';

  fin = fopen(rawfilename, "rb");
  gain = 1000;
  frames = nbits/(Nc*Nb);

  prev_rx_symbols = ones(Nc+1,1);
  foff_phase = 1;

  % BER stats

  total_bit_errors = 0;
  total_bits = 0;
  bit_errors_log = [];
  sync_log = [];
  test_frame_sync_log = [];
  test_frame_sync_state = 0;
  error_pattern_log = [];

  % SNR states

  sig_est = zeros(Nc+1,1);
  noise_est = zeros(Nc+1,1);

  % logs of various states for plotting

  rx_symbols_log = [];
  rx_timing_log = [];
  foff_log = [];
  rx_fdm_log = [];
  snr_est_log = [];

  % misc states

  nin = M; % timing correction for sample rate differences
  foff = 0;
  track_log = [];
  track = 0;
  fest_state = 0;

  % spectrum states

  Nspec=1024;
  spec_mem=zeros(1,Nspec);
  SdB = zeros(1,Nspec);

  % Main loop ----------------------------------------------------

  for f=1:frames
    
    % obtain nin samples of the test input signal
    
    for i=1:nin
      rx_fdm(i) = fread(fin, 1, "short")/gain;
    end
    
    rx_fdm_log = [rx_fdm_log rx_fdm(1:nin)];

    % update spectrum

    l=length(rx_fdm);
    spec_mem(1:Nspec-l) = spec_mem(l+1:Nspec);
    spec_mem(Nspec-l+1:Nspec) = rx_fdm;
    S=fft(spec_mem.*hanning(Nspec)',Nspec);
    SdB = 0.9*SdB + 0.1*20*log10(abs(S));

    % frequency offset estimation and correction

    [pilot prev_pilot pilot_lut_index prev_pilot_lut_index] = get_pilot(pilot_lut_index, prev_pilot_lut_index, nin);
    [foff_coarse S1 S2] = rx_est_freq_offset(rx_fdm, pilot, prev_pilot, nin);
    
    if track == 0
      foff  = foff_coarse;
    end
    foff_log = [ foff_log foff ];
    foff_rect = exp(j*2*pi*foff/Fs);

    for i=1:nin
      foff_phase *= foff_rect';
      rx_fdm(i) = rx_fdm(i)*foff_phase;
    end

    % baseband processing

    rx_baseband = fdm_downconvert(rx_fdm, nin);
    rx_filt = rx_filter(rx_baseband, nin);

    [rx_symbols rx_timing] = rx_est_timing(rx_filt, rx_baseband, nin);
    
    rx_timing_log = [rx_timing_log rx_timing];
    nin = M;
    if rx_timing > 2*M/P
       nin += M/P;
    end
    if rx_timing < 0;
       nin -= M/P;
    end

    if strcmp(modulation,'dqpsk')
      rx_symbols_log = [rx_symbols_log rx_symbols.*conj(prev_rx_symbols./abs(prev_rx_symbols))*exp(j*pi/4)];
    else
      rx_symbols_log = [rx_symbols_log rx_symbols];
    endif
    [rx_bits sync f_err pd] = qpsk_to_bits(prev_rx_symbols, rx_symbols, modulation);
    [sig_est noise_est] = snr_update(sig_est, noise_est, pd);
    snr_est = calc_snr(sig_est, noise_est);
    snr_est_log = [snr_est_log snr_est];
    foff -= 0.5*f_err;
    prev_rx_symbols = rx_symbols;
    sync_log = [sync_log sync];

    % freq est state machine

    [track fest_state] = freq_state(sync, fest_state);
    track_log = [track_log track];

    % count bit errors if we find a test frame

    [test_frame_sync bit_errors error_pattern] = put_test_bits(test_bits, rx_bits);
    if (test_frame_sync == 1)
      total_bit_errors = total_bit_errors + bit_errors;
      total_bits = total_bits + Ntest_bits;
      bit_errors_log = [bit_errors_log bit_errors/Ntest_bits];
    else
      bit_errors_log = [bit_errors_log 0];
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
        else
          error_pattern_log = [error_pattern_log error_pattern];
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
  Fcentre
  printf("%d bits  %d errors  BER: %1.4f\n",total_bits, total_bit_errors, ber);

  % ---------------------------------------------------------------------
  % Plots
  % ---------------------------------------------------------------------

  xt = (1:frames)/Rs;
  secs = frames/Rs;

  figure(1)
  clf;
  [n m] = size(rx_symbols_log);
  plot(real(rx_symbols_log(1:Nc+1,15:m)),imag(rx_symbols_log(1:Nc+1,15:m)),'+')
  axis([-2 2 -2 2]);
  title('Scatter Diagram');

  figure(2)
  clf;
  subplot(211)
  plot(xt, rx_timing_log)
  title('timing offset (samples)');
  subplot(212)
  plot(xt, foff_log, '-;freq offset;')
  hold on;
  plot(xt, track_log*75, 'r;course-fine;');
  hold off;
  title('Freq offset (Hz)');
  grid

  figure(3)
  clf;
  spec(rx_fdm_log,8000);

  figure(4)
  clf;
  subplot(311)
  stem(xt, sync_log)
  axis([0 secs 0 1.5]);
  title('BPSK Sync')
  subplot(312)
  stem(xt, bit_errors_log);
  title('Bit Errors for test frames')
  subplot(313)
  plot(xt, test_frame_sync_log);
  axis([0 secs 0 1.5]);
  title('Test Frame Sync')

  figure(5)
  clf;
  subplot(211);
  plot(xt, snr_est_log);
  title('SNR Estimates')
  subplot(212)
  snrdB_pc = 20*log10(sig_est(1:Nc+1)) - 20*log10(noise_est(1:Nc+1));
  bar(snrdB_pc(1:Nc) - mean(snrdB_pc(1:Nc)))
  axis([0 Nc+1 -3 3]);

  figure(6)
  clf;
  hold on;
  lep = length(error_pattern_log);
  for p=1:Nc
    plot(p + 0.25*error_pattern_log((p-1)*2+1:Nc*Nb:lep));
    plot(0.30 + p + 0.25*error_pattern_log(p*2:Nc*Nb:lep),'r')
  end
  hold off;
  axis([1 lep/(Nc*Nb) 0 15])

  figure(7)
  clf;
  subplot(211)
  [a b] = size(rx_fdm_log);
  xt1 = (1:b)/Fs;
  plot(xt1, rx_fdm_log);
  title('Rx FDM Signal');
  subplot(212)
  plot((0:Nspec/2-1)*Fs/Nspec, SdB(1:Nspec/2) - 20*log10(Nspec/2))
  axis([0 Fs/2 -40 0])
  grid
  title('FDM Rx Spectrum');

endfunction
