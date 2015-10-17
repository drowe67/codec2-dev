% legacyfsk.m
% David Rowe October 2015
%
% Attempt at FSK waveform that can pass through legacy FM radios
% but still be optimally demodulated by SDRs. It doesn't have to be
% optimally demodulated by legacy radios.  Trick is getting it to
% pass through 300-3000Hz audio filters in leagcy radios.
%
% [X] code up modulator
%     [X] manchester two bit symbols
%     [X] plot spectrum
% [ ] demodulate
% [ ] measure BER compared to ideal coherent FSK

1;

fm; % analog FM library


function states = legacyfsk_init()
  Fs = states.Fs = 96000;
  Rs = states.Rs = 4800;          
  Ts = states.Ts = Fs/Rs;
  nbits = states.nbits = 100;                  % number of payload data symbols/frame
  nbits2 = states.nbits2 = states.nbits*2;     % number of bits/frame over channel after manchester encoding
endfunction


% test modulator function

function tx = legacyfsk_mod(states, tx_bits)
    tx = zeros(states.Ts*length(tx_bits),1);
    tx_phase = 0;
    Ts = states.Ts;
    Fs = states.Fs;
    Rs = states.Rs;
    f1 = 24E3-Rs/2; f2 = 24E3+Rs/2;

    for i=1:length(tx_bits)
      if tx_bits(i) == 0
        tx_phase_vec = tx_phase + (1:Ts)*2*pi*f1/Fs;
      else
        tx_phase_vec = tx_phase + (1:Ts)*2*pi*f2/Fs;
      end
      tx((i-1)*Ts+1:i*Ts) = 2.0*cos(tx_phase_vec);
      tx_phase = tx_phase_vec(Ts) - floor(tx_phase_vec(Ts)/(2*pi))*2*pi;
    end
endfunction


function run_sim

  frames = 10;
  EbNodB = 12.5;
  test_frame_mode = 1;

  % init fsk modem

  more off
  rand('state',1); 
  randn('state',1);
  states = legacyfsk_init();
  Fs = states.Fs;
  nbits = states.nbits;
  nbits2 = states.nbits2;
  Ts = states.Ts;

  % init analog FM modem

  fm_states.Fs = Fs;  
  fm_max = fm_states.fm_max = 3E3;
  fd = fm_states.fd = 5E3;
  fm_states.fc = 24E3;

  fm_states.pre_emp = 0;
  fm_states.de_emp  = 1;
  fm_states.Ts = 1;
  fm_states.output_filter = 1;
  fm_states = analog_fm_init(fm_states);
  [b, a] = cheby1(4, 1, 300/Fs, 'high');   % 300Hz HPF to simulate FM radios

  rx_bits_buf = zeros(1,2*nbits);
  Terrs = Tbits = 0;
  state = 0;
  nerr_log = [];

  EbNo = 10^(EbNodB/10);
  variance = states.Fs/((states.Rs/2)*EbNo);  % actual bit rate is Rs/2

  % manchester code templates (matched filter coefficients)

  manchester_one = ones(1,2*Ts);
  manchester_zero = ones(1,2*Ts);
  manchester_one(Ts+1:2*Ts) = manchester_zero(1:Ts) = -1;

  if test_frame_mode == 1
    % test frame of bits, which we repeat for convenience when BER testing
    test_frame = round(rand(1, states.nbits));
    tx_bits = [];
    for i=1:frames+1
      tx_bits = [tx_bits test_frame];
    end
  end
  if test_frame_mode == 2
    % random bits, just to make sure sync algs work on random data
    tx_bits = round(rand(1, states.nbits*(frames+1)));
  end
  if test_frame_mode == 3
    % ...10101... sequence
    tx_bits = zeros(1, states.nbits*(frames+1));
    tx_bits(1:2:length(tx_bits)) = 1;
    %tx_bits(10:length(tx_bits)) = 1;
  end

  % Manchester encode, halving the payload bit rate, and removing DC
  % term in baseband signal, which makes it a bit more friendly to 
  % old-school legacy FM radios.

  tx_bits_mapped = zeros(1,length(tx_bits)*2);
  j = 1;
  for i=1:2:length(tx_bits_mapped)
    if tx_bits(j)
      tx_bits_mapped(i) = 1;
      tx_bits_mapped(i+1) = 0;
    else
      tx_bits_mapped(i) = 0;
      tx_bits_mapped(i+1) = 1;
    end
    j++;
  end

  % use ideal FSK modulator (note: need to try using analog FM modulator)

  tx = legacyfsk_mod(states, tx_bits_mapped);
  noise = sqrt(variance)*randn(length(tx),1);
  rx    = tx + noise;

  % use analog FM demodulator

  [rx_out rx_bb] = analog_fm_demod(fm_states, rx');
  rx_out_hp = filter(b,a,rx_out);

  % filter using manchester code templates, aka integration or matched filter

  rx_filt_one = filter(manchester_one,1,rx_out_hp);
  rx_filt_zero = filter(manchester_zero,1,rx_out_hp);
  
  rx_timing = 5;
  rx_filt_one_dec = rx_filt_one(rx_timing:2*Ts:length(rx_filt_one));
  rx_filt_zero_dec = rx_filt_zero(rx_timing:2*Ts:length(rx_filt_zero));

  st = 1;
  for f=1:frames

    % extract nin bits

    nin = nbits;
    en = st + nin - 1;
    rx_bits = rx_filt_one_dec(st:en) < rx_filt_zero_dec(st:en);
    st += nin;

    rx_bits_buf(1:nbits) = rx_bits_buf(nbits+1:2*nbits);
    rx_bits_buf(nbits+1:2*nbits) = rx_bits;

    % frame sync based on min BER

    if test_frame_mode == 1
      nerrs_min = nbits;
      next_state = state;
      if state == 0
        for i=1:nbits
          error_positions = xor(rx_bits_buf(i:nbits+i-1), test_frame);
          nerrs = sum(error_positions);
          % printf("i: %d nerrs: %d nerrs_min: %d \n", i, nerrs, nerrs_min);
          if nerrs < nerrs_min
            nerrs_min = nerrs;
            coarse_offset = i;
          end
        end
        if nerrs_min < 3
          next_state = 1;
          %printf("%d %d\n", coarse_offset, nerrs_min);
        end
      end

      if state == 1  
        error_positions = xor(rx_bits_buf(coarse_offset:coarse_offset+nbits-1), test_frame);
        nerrs = sum(error_positions);
        Terrs += nerrs;
        Tbits += nbits;
        nerr_log = [nerr_log nerrs];
      end

      state = next_state;

    end 
  end

  if test_frame_mode == 1
    printf("frames: %d Tbits: %d Terrs: %d BER %4.3f\n", frames, Tbits, Terrs, Terrs/Tbits);
  end

  % Bunch O'plots --------------------------------------------------------------

  Tx=fft(tx, Fs);
  TxdB = 20*log10(abs(Tx(1:Fs/2)));
  figure(1)
  clf;
  plot(TxdB)
  axis([1 Fs/2 (max(TxdB)-100) max(TxdB)])
  title('Tx Spectrum');

  figure(2)
  clf
  subplot(211)
  st = 1; en=20;
  plot(rx_out(st:en*states.Ts*2));
  title('After Analog FM demod');
  subplot(212)
  plot(rx_out_hp(st:en*states.Ts*2));
  title('After 300Hz HPF');

  figure(3);
  clf;
  subplot(211)
  h = freqz(b,a,Fs);
  plot(20*log10(abs(h(1:4000))))
  title('300Hz HPF Response');
  subplot(212)
  h = fft(rx_out, Fs);
  plot(20*log10(abs(h(1:4000))))
  title('Spectrum of baseband modem signal after analog FM demod');

  figure(4);
  clf;
  subplot(211)
  plot(rx_filt_one(st:en*states.Ts*2) - rx_filt_zero(st:en*states.Ts*2),'+')
  title('Difference after matched filtering');
  subplot(212)
  plot(rx_filt_one_dec(st:en),'+');
  hold on;
  plot(rx_filt_zero_dec(st:en),'g+')
  hold off;
  title('Both channels after sampling at ideal timing instant')
  
  hold off;
end

run_sim;
