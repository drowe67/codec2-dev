% impulse.m
% experimenting with impulse spread over time

N=256;
Fs=8000;
f0=100;

function [s_real phi env] = synth(disperse=0, filt=0, f0, Fs, N)
  Wo=2*pi*f0/Fs;
  L=floor(pi/Wo);
  n=0:N-1;
  n0=[0 1 2 4];
  s = zeros(1,N);
  phi=zeros(1,L);
  for k=1:4
    st = 1 + floor((k-1)*L/4);
    en = floor(k*L/4);
    for m=st:en
      if disperse; phi(m) = -Wo*m*n0(k); end
      s += exp(j*(Wo*m*n + phi(m)));
    end
  end
  beta = 0.9;
  s_real = real(s);
  env = abs(s);
  if filt; s_real = filter(1,[1 -2*beta*cos(pi/4) beta*beta], s_real); end
endfunction

[s_p0 phi_p0 env_p0] = synth(disperse=0, filt=0, f0, Fs, N);
save_raw("imp_p0.raw",1000*s_p0);
[s_disp phi_disp env_disp] = synth(disperse=1, filt=0, f0, Fs, N);
save_raw("imp_disp.raw",1000*s_filt);

figure(1); clf;
subplot (211); plot(s_p0); hold on; plot(env_p0,'g'); hold off; %axis([1 200 -100 100]);
subplot (212); plot(s_disp); hold on; plot(env_disp,'g'); hold off; %axis([1 200 -100 100]);
figure(2); clf;
subplot(211); plot(phi_disp,"o-");
Wo=2*pi*f0/Fs;
L=floor(pi/Wo);
phase_delay = phi_disp./(Wo*(1:L)*Fs/(2*pi));
group_delay = (phi_disp(2:L) - phi_disp(1:L-1))/(Wo*Fs/(2*pi));
subplot(212); plot(group_delay,"o-");
