% legacyfsk.m
% David Rowe October 2015
%
% Attempt at FSK waveform that can pass through legacy FM radios
% but still be optimally demodulated by SDRs. It doesn't have to be
% optimally demodulated by legacy radios.  Trick is getting it to
% pass through 300-3000Hz audio filters in leagcy radios.
%
% [ ] code up modulator
%     [ ] map to two bit symbols
%     [ ] plot spectrum
% [ ] can we just get away with a scrambler?

1;

function states = legacyfsk_init()
  Fs = states.Fs = 96000;
  Rs = states.Rs = 4800;          
  Ts = states.Ts = Fs/Rs;
endfunction


% test modulator function

function tx = legacyfsk_mod(states, tx_bits)
    tx = zeros(states.Ts*length(tx_bits),1);
    tx_phase = 0;
    Ts = states.Ts;
    Fs = states.Fs;
    Rs = states.Rs
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


fm;

states = legacyfsk_init();
nbits = states.Rs;
Fs= states.Fs;
test_mode = 1;
if test_mode == 1
  tx_bits = round(rand(1, nbits/2));
else
  % ...10101... sequence
  tx_bits = zeros(1, nbits/2);
  tx_bits(1:2:length(tx_bits)) = 1;
end

tx_bits_mapped = zeros(1,nbits);
j = 1;
for i=1:2:nbits
  if tx_bits(j)
    tx_bits_mapped(i) = 1;
    tx_bits_mapped(i+1) = 0;
  else
    tx_bits_mapped(i) = 0;
    tx_bits_mapped(i+1) = 1;
  end
  j++;
end

tx = legacyfsk_mod(states, tx_bits_mapped);
Tx=fft(tx);
TxdB = 20*log10(abs(Tx(1:Fs/2)));
figure(1)
clf;
plot(TxdB)
axis([1 Fs/2 (max(TxdB)-100) max(TxdB)])

fm_states.Fs = Fs;  
fm_max = fm_states.fm_max = 3E3;
fd = fm_states.fd = 5E3;
fm_states.fc = 24E3;

fm_states.pre_emp = 0;
fm_states.de_emp  = 1;
fm_states.Ts = 1;
fm_states.output_filter = 0;
fm_states = analog_fm_init(fm_states);

[rx_out rx_bb] = analog_fm_demod(fm_states, tx');
[b, a] = cheby1(4, 1, 300/Fs, 'high');
rx_out_hp = filter(b,a,rx_out);

figure(2)
clf
subplot(211)
st = 1; en=40*states.Ts;
plot(rx_out(st:en));
subplot(212)
plot(rx_out_hp(st:en) > 0);
axis([1 en -0.5 1.5])

figure(3);
clf;
subplot(211)
h = freqz(b,a,Fs);
plot(20*log10(abs(h(1:4000))))
subplot(212)
h = fft(rx_out);
plot(20*log10(abs(h(1:4000))))
