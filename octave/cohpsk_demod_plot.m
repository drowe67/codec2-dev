% cohpsk_demod_plot.m
% David Rowe May 2015
%
% Plot Octave outputs from cohpsk_demod, c2dec, to visualise whats going on
% when errors hit the system

% $ ./c2enc 700 ../../raw/ve9qrp_10s.raw - | ./cohpsk_mod - - | ./cohpsk_ch - - -60 50 1 1 | ./cohpsk_demod - - cohpsk_demod.txt | ./c2dec 700 - - --dump ve9qrp | play -t raw -r 8000 -s -2 - -q

% ./c2enc 700 ../../raw/ve9qrp_10s.raw - | ./cohpsk_mod - - | ./cohpsk_ch - - -30 50 1 1 | ./cohpsk_demod - - cohpsk_demod.txt | ./c2dec 700 - - --dump ve9qrp_snr3 | play -t raw -r 8000 -s -2 - -q

graphics_toolkit ("gnuplot");

Nc=7; Nd=2; Ns=6;

load ../build_linux/src/cohpsk_demod.txt
load ../build_linux/src/cohpsk_put_test_bits.txt
load ../build_linux/src/ve9qrp_lsp_.txt
load ../build_linux/src/ve9qrp_snr3_lsp_.txt
load ../build_linux/src/ve9qrp_ak_.txt
load ../build_linux/src/ve9qrp_snr3_ak_.txt
load ../build_linux/src/ve9qrp_model.txt
load ../build_linux/src/ve9qrp_snr3_model.txt
load ../build_linux/src/ve9qrp_snr3_softdec.txt

Ncf = 50;     % number of codec frames to plot
Nmf = Ncf/2;  % number of modem frames to plot
Nms = Nmf*Ns; % number of modem symbols to plot

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
plot(rx_phi_log_c(1:Nms,:))
title('phase')
axis([1 Nms -pi pi])
subplot(212)
plot(rx_amp_log_c(1:Nms,:))
title('amplitude')
axis([1 Nms 0 1])

figure(3)
subplot(211)
plot(rx_timing_log_c)
title('rx timing');
subplot(212)
stem(ratio_log_c)
title('Sync ratio');

figure(4);
clf;
plot(nerr_log_c(1:Ncf));
title('Bit Errors');
xlabel('Codec Frame')

figure(5);
clf;
plot(error_positions_hist_c);
title('Error Position Histogram');

figure(6)
y = 1:Nms;
x = 1:Nc*Nd;
z = 20*log10(rx_amp_log_c(1:Nms,:));
mesh(x,y,z);
grid
title('Channel Amplitude dB');
a = min(min(z));
b = max(max(z));
axis([1 Nc*Nd 1 Nms a b])

% work out alignment, as they sync at different times

min_e = 1E6;
for i=1:10
  l1 = length(ve9qrp_lsp_);
  l2 = length(ve9qrp_snr3_lsp_);
  st = i; en = min(l1+i-1,l2);
  d = ve9qrp_lsp_(st:en, 1:6) - ve9qrp_snr3_lsp_(1:en-st+1, 1:6);
  e = sum(sum(abs(d)));
  if e < min_e
    min_e = e;
    min_i = i;
  end
end
printf("time offset between clean and 3dB is %d codec frames\n", min_i);

% LSP trajectories

figure(7)
clf
st = min_i; en = 50;
plot(ve9qrp_snr3_lsp_(1:en-st+1, 1:6),'r--')
hold on
plot(ve9qrp_lsp_(st:en, 1:6),'g')
hold off

% Spectral distortion of LPCs

figure(8)
clf;
f1=1./fft(ve9qrp_ak_(st:en,:)',128);
f2=1./fft(ve9qrp_snr3_ak_(1:en-st+1,:)',128);
%d = (20*log10(abs(f1)) - 20*log10(abs(f2)));
d = 20*log10(abs(f2));
sdsq = mean(d.^2);
plot(sdsq)
title('spectral distortion clean and channel SNR=3dB')

figure(9)
clf;
y = 1:en-st+1;
x = 1:40;
%mesh(y,x,-20*log10(abs(f2(1:40,:))));
mesh(y,x,d(x,:));
grid
title('Synthesis filter difference between clean and channel SNR=3dB');
xlabel('Time (codec frames)')
ylabel('Frequency 0 to 2500Hz');
zlabel('Difference (dB)');

% map soft decn information to LSPs

mel1 = ve9qrp_snr3_softdec(:,10:12);
mel2 = ve9qrp_snr3_softdec(:,13:14);
mel3 = ve9qrp_snr3_softdec(:,15:18);
mel4 = ve9qrp_snr3_softdec(:,19:21);
mel5 = ve9qrp_snr3_softdec(:,22:24);
mel6 = ve9qrp_snr3_softdec(:,25:26);
softdec_mel = [sum(mel1'.^2); sum(mel2'.^2); sum(mel3'.^2); sum(mel4'.^2); sum(mel5'.^2); sum(mel6'.^2)];

figure(10)
clf;
y = 1:en-st+1;
x = 1:6;
%mesh(y,x,-20*log10(abs(f2(1:40,:))));
mesh(y, x, softdec_mel(:,y));
grid
xlabel('Codec frame')
ylabel('LSP')
zlabel('Power')
%axis([1 (en-st+1) 1 6 -10 5])

% plot symbol energy against SD

figure(11)
semilogx(mean(softdec_mel(:,1:en-st+1)), sdsq,'+')
grid
