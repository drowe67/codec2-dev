% impulse.m
% experimenting with impulse spread over time

N=8000;
Fs=8000;
f0=80;

function [s_sum phi s] = synth(disperse=0, filt=0, f0, Fs, N)
  Wo=2*pi*f0/Fs;
  L=floor(pi/Wo);
  n=0:N-1;
  n0=[0 1 2 4];
  s = zeros(4,N);
  phi=zeros(1,L);
  for k=1:4
    st = 1 + floor((k-1)*L/4);
    en = floor(k*L/4);
    for m=st:en
      if disperse; phi(m) = -Wo*m*n0(k); end
      s(k,:) += cos(Wo*m*n + phi(m));
    end
  end
  s_sum = sum(s);
  beta = 0.9;
  if filt; s_sum = filter(1,[1 -2*beta*cos(pi/8) beta*beta],s_sum); end
endfunction

s_p0 = synth(disperse=0, filt=0, f0, Fs, N);
save_raw("imp_p0.raw",1000*s_p0);
[s_filt phi] = synth(disperse=1, filt=0, f0, Fs, N);
save_raw("imp_disp.raw",1000*s_filt);

figure(1); clf;
subplot (211); plot(s_p0); axis([1 200 -50 50]);
subplot (212); plot(s_filt); axis([1 200 -50 50]);
figure(2); clf;
subplot(211); plot(phi,"o-");
Wo=2*pi*f0/Fs;
L=floor(pi/Wo);
phase_delay = phi./(Wo*(1:L)*Fs/(2*pi));
group_delay = (phi(2:L)-phi(1:L-1))/(Wo*Fs/(2*pi));
subplot(212); plot(group_delay,"o-");
