% cohpsk_demod.m
% David Rowe May 2015
%
% Plot Octave outputs from cohpsk_demod_plot

Nc=7; Nd=2;

load ../build_linux/src/cohpsk_demod.txt
load ../build_linux/src/cohpsk_put_test_bits.txt
  
figure(1)
clf;

% plot combined signals to show diversity gains

combined = rx_symb_log_c(:,1:Nc);
for d=2:Nd
  combined += rx_symb_log_c(:, (d-1)*Nc+1:d*Nc);
end
plot(combined*exp(j*pi/4)/sqrt(Nd),'+')
title('Scatter');
axis([-2 2 -2 2])

figure(2)
clf;
subplot(211)
plot(rx_phi_log_c)
title('phase')
subplot(212)
plot(rx_amp_log_c)
title('amplitide')

figure(3)
subplot(211)
plot(rx_timing_log_c)
title('rx timing');
subplot(212)
stem(ratio_log_c)
title('Sync ratio');

figure(4);
clf;
plot(nerr_log_c);
title('Bit Errors');

figure(5);
clf;
plot(error_positions_hist_c);
title('Error Position Histogram');
