% impulse.m
% experimenting with impulse spread over time

N=256;
Fs=8000;
f0=80;

function t0 = find_t0(dt0,df,Tspread,f)
  t0 = mod( (dt0/df)*f, Tspread);
endfunction

function test_t0
  figure(3); clf;
  f = 100:100:1200;
  t0 = [];
  for i=1:length(f)
    af = f(i);
    t0 = [t0 find_t0(0.5E-3, 100, 2E-3, af)];
  end
  plot(f,t0);
endfunction

function [s_real phi env] = synth(disperse=0, filt=0, f0, Fs, N)
  Wo=2*pi*f0/Fs;
  L=floor(pi/Wo);
  n=0:N-1;
  s = zeros(1,N);
  phi = zeros(1,L);
  if disperse
    for m=1:L
      f = m*Wo*Fs/(2*pi);
      t0m = find_t0(0.1E-3, 100, 1E-3, f);
      phi(m) = -Wo*m*t0m*Fs;
    end
  end
  for m=1:L
    s += exp(j*(Wo*m*n + phi(m)));
  end
  s_real = real(s);
  env = abs(s);
  beta = 0.95;
  if filt
    s_real = filter(1,[1 -2*beta*cos(pi/6) beta*beta], s_real);
  end
endfunction

[s_p0 phi_p0 env_p0] = synth(disperse=0, filt=1, f0, Fs, N);
save_raw("imp_p0.raw",1000*s_p0);
[s_disp phi_disp env_disp] = synth(disperse=1, filt=1, f0, Fs, N);
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
subplot(212); plot(phase_delay,"o-");

test_t0
