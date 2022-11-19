% tdetphase.m
% David Rowe August 2017
%
% Testing Hilbert Transform recover of phase from magnitude spectra

newamp_700c;
Fs = 8000;

w = 2*pi*500/Fs; gamma = 0.97;
ak = [1 -2*gamma*cos(w) gamma*gamma];
Nfft = 512;

% Test 1 - compare phase from freqz for 2nd order system (all pole filter)
%        - uses internal test of determine_phase()

h = freqz(1,ak,Nfft/2);

% note dummy_model not used, as determine_phase() is used in test mode

L = 20; Wo = pi/(L+1);
dummy_model = [Wo L ones(1,L)];
phase = determine_phase(dummy_model, 1, Nfft, ak);

fg = 1;
figure(fg++); clf;
subplot(211); plot(20*log10(abs(h))); ylabel('Amplitude');
title('Test 1 - phase from high sample rate');
subplot(212); plot(angle(h)); hold on; plot(phase(1:Nfft/2),'g+'); hold off; ylabel('Phase');

test1_error = sum((angle(h) - phase(1:Nfft/2)(1:256)').^2);
if (test1_error < 1E-2), printf('Test 1 PASS\n'); else printf('Test 1 FAIL\n'); end

% Test 2 - Use phase samples at harmonics

F0 = 100; Wo = 2*pi*F0/Fs; L = floor(pi/Wo);
Am = zeros(1,L);
for m=1:L
  b = round(m*Wo*Nfft/(2*pi));
  Am(m) = abs(h(b));
end
AmdB = 20*log10(Am);
model = [Wo L Am];
[phase Gdbfk s] = determine_phase(model, 1, Nfft);

fftx = (0:Nfft/2-1)*(Fs/Nfft);
harmx = (1:L)*Wo*Fs/(2*pi);

figure(fg++); clf;
subplot(211); plot(fftx, Gdbfk(1:Nfft/2)); ylabel('Gdbfk');
title('Test 2 - Gdbfk and s from determine\_phase()');
subplot(212); plot(s(1:Nfft/2)); ylabel('s');

figure(fg++); clf;
subplot(211); plot(fftx, 20*log10(abs(h(1:Nfft/2))));
              hold on; plot(harmx, AmdB, 'g+;orig;','markersize',10,'linewidth',2);
              plot(fftx, Gdbfk(1:Nfft/2), 'r;synth;'); hold off;
title('Test 2 - phase from F0=100Hz harmonics');
subplot(212); plot(fftx, angle(h), 'g;orig;'); hold on; plot(fftx, phase(1:Nfft/2),'r;synth;'); hold off;

# note this fails - possibly due to aliasing or para interpolation.  Shows issues with synth phase from
# discrete spectra

test2_error = sum((angle(h(1:Nfft/2)) - phase(1:Nfft/2)(1:256)').^2);
if (test2_error < 1E-2), printf('Test 2 PASS\n'); else printf('Test 2 FAIL (as expected)\n'); end

#{
  Test 3 - optional demo using real harmonic amplitudes

  Create inut files with
  $ cd ~/codec2/buildlinux
  $ ./src/c2sim ../raw/hts1a.raw --phase0 --dump hts1a
#}

if ((exist("../build_linux/hts1a_model.txt") == 2) && (exist("../build_linux/hts1a_phase.txt") == 2))
  model = load("../build_linux/hts1a_model.txt");
  phase_orig = load("../build_linux/hts1a_phase.txt");

  f = 42;
  Wo = model(f,1); L = model(f,2); Am = model(f,3:L+2); AmdB = 20*log10(Am);
  [phase Gdbfk s] = determine_phase(model, f, Nfft);

  fftx = (1:Nfft/2)*(Fs/Nfft);
  harmx = (1:L)*Wo*Fs/(2*pi);

  figure(fg++); clf;
  subplot(211); plot(fftx, Gdbfk(1:Nfft/2));
  subplot(212); plot(s(1:Nfft/2))

  figure(fg++); clf;
  subplot(211); plot(harmx, AmdB, 'g+;AmdB;','markersize',10);
                hold on; plot(fftx, Gdbfk(1:Nfft/2), 'r;Gdbfk;'); hold off;
                axis([0 Fs/2 0 80]);
  subplot(212); plot(fftx, phase(1:Nfft/2),'g;HT phase;');

  % synthesise using phases

  N = 320;
  s = s_phase = zeros(1,N);
  for m=1:L
    s = s + Am(m)*cos(m*Wo*(1:N) + phase_orig(f,m));
    b = round(m*Wo*Nfft/(2*pi));
    s_phase = s_phase + Am(m)*cos(m*Wo*(1:N) + phase(b));
  end
  figure(fg++); clf;
  subplot(211); plot(s); title('Speech from Orig Phase');
  subplot(212); plot(s_phase,'g'); title('Speech From HT Phase');
end
