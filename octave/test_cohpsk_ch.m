% test_cohpsk_ch.m
% David Rowe May 2015
%
% Plot outputs from test_coh_psk_ch.c

Nc=7; Nd=2;

load ../build_linux/src/test_cohpsk_ch_out.txt
  
figure(3)
clf;

% plot combined signals to show diversity gains

combined = rx_symb_log_c(:,1:Nc);
for d=2:Nd
  combined += rx_symb_log_c(:, (d-1)*Nc+1:d*Nc);
end
plot(combined*exp(j*pi/4)/sqrt(Nd),'+')
title('Scatter');
axis([-2 2 -2 2])
