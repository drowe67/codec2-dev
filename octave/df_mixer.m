% df_mixer.m
%
% David Rowe October 2015
%
% Experimental direction finding using a mixer to FDM signals from two
% antennas.  This simulation tests estimation of phase and freq of 
% mixer LO.
%
% [X] built in hackrf_dc
% [X] test mode
% [X] compass needle
% [X] real time operation (or close to it)

fm;

% Generate a FM modulated signal for testing ----------------------------------

function tx = fm_mod(Fs, fmod)
  fm_states.Fs = Fs;  
  fm_states.fm_max = 3E3;
  fm_states.fd = 5E3;
  fm_states.fc = 48E3;
  fm_states.pre_emp = 0;
  fm_states.de_emp  = 0;
  fm_states.Ts = 1;
  fm_states.output_filter = 1;
  fm_states = analog_fm_init(fm_states);

  nsam = fm_states.Fs/2;
  t = 0:(nsam-1);
  wmod = 2*pi*fmod/fm_states.Fs;
  mod = sin(wmod*t);
  tx = analog_fm_mod(fm_states, mod);
endfunction


% Check Filters -------------------------------------------------------------

% Run some test signals through all three filters.  Real signals used
% for convenience, also works with complex.  If another set of filter
% is used you need to prove net reponse is a delay.  Fig 2
% subplot(211) must be a straight line.

function check_filters(bcentre, Fs, f, bw, w, wlo, t)
  x = cos(w*t) + cos((w-wlo)*t) + cos((w-wlo)*t); 
  ycenter = filter(bcentre,1,x);
  ylow    = exp(-j*wlo*t) .* filter(bcentre, 1, x .* exp(j*wlo*t));
  yhigh   = exp(j*wlo*t) .* filter(bcentre, 1, x .* exp(-j*wlo*t));
  y = ycenter + ylow + yhigh;

  figure(1)
  clf
  h=freqz(bcentre,1,Fs/2);
  plot(20*log10(abs(h)))
  hold on;
  plot([f-bw/2 f-bw/2 f+bw/2 f+bw/2],[-60 -3 -3 -60],'r');
  hold off;
  axis([f-bw f+bw -60 0])
  title('Centre BPF');
  grid;

  figure(2)
  clf;
  Ns = 200;
  Np = 200;

  % y should be a delayed version of x

  subplot(211)
  plot(x(Ns:Ns+Np-1))
  hold on;
  plot(y(Ns:Ns+Np-1),'g')
  hold off;
  title('Filter Check');

  % should see a straight line on this plot

  subplot(212)
  plot(x(Ns:Ns+Np-1), y(Ns:Ns+Np-1))
end


% Constants -------------------------------------------------------------------

randn('state',1);
more off;
format short

f     = 48E3;          % signal frequency
flo   = 32.768E3;      % LO frequency
m     = 10^(-0/20);    % modulation index of our mixer, how far each 
                       % sideband is down wrt carrier (0.1 == 20dB)
                       % note: not FM modulation index
bw    = 16E3;          % signal bandwidth

phi   = 0.5;           % phase difference between antennas - 
                       % ->> what we are trying to estimate

theta = 0;             % phase term constant to ant1 and ant2
alpha = 0;             % LO phase

Fs    = 200E3;         % sample rate
Fshrf = 10E6;          % HackRF sample rate
CNodB = 100;           % C/No of carrier, which dominates power
                       % Around C/No = 50dB is min for FM (12dB-ish SINAD)
fmod  = 0;             % Audio modulation on FM signal (0 for carrier)
N     = Fs/2;          % processing block size

% BW is Fs Hz.  No=(total noise power Nt)/Fs. N=noise power=variance of noise source
% C = carrier power = 1
% CNo = C/No = 1/(Nt/Fs), therefore var = Nt = Fs/CNo

CNo = 10^(CNodB/10);
var = Fs/CNo;

w   = 2*pi*f/Fs;
wlo = 2*pi*flo/Fs;
wbw = 2*pi*bw/Fs;

% FIR filter BPF for centre signal.  Use a 15kHz cut off to catch
% energy from FM signals.  Important: we need identical phase response
% for each BPF.  I couldn't work out how to design these filter using
% the built-in Octave functions.  So I use a bunch of freq shifts at
% run time instead, so we use the centre BPF for all three filters.

bcentre  = fir2(200, [0, f-bw/2, f, f+bw/2, Fs/2]/(Fs/2), [0 0 1 0 0]);

%check_filters(bcentre, Fs, f, bw, w, wlo, t); xx

% Main -------------------------------------------------------------------------

mode = "sample";
finished = 0;
while !finished

  if strcmp(mode, "simulate")
    % Simulate rx signal ---------------------------------------------------------

    t=0:N-1;

    % tx signal

    tx = fm_mod(Fs, fmod); 

    % simulate rx signals at antenna 1 and 2 by adding noise.  Note signal
    % at antenna 2 is phase shifted phi due to path difference. Both signals
    % experience a time of flight path different theta.

    n = sqrt(var/2)*randn(1,N) + j*sqrt(var/2)*randn(1,N);
    ant1 = tx*exp(j*(theta+phi)) + n;
    ant2 = tx*exp(j*(theta)) + n;

    % ant2 passed through mixer with local osc frequency wlo and phase alpha

    ant2_mix = 2.0*m*(ant2 .* cos(t*wlo + alpha));

    % The SDR receives the sum of the two signals

    rx = ant1 + ant2_mix;
  end

  if strcmp(mode,"file") || strcmp(mode,"sample")
    % Off air signal from HackRF stored in file ----------------------------------------------

    if strcmp(mode,"sample")
      [status output] = system("hackrf_transfer -r df1.iq -f 439000000  -n 10000000 -l 20 -g 40");
    end
    s1 = load_hackrf("df1.iq");
    rx = downsample(rot90(s1), Fshrf/Fs)/127;
    t = 1:length(rx);  
  end

  Rx = (1/N)*fft(rx,N);

  % BPF each signal

  sam1 = filter(bcentre,1,rx);
  sam2 = exp( j*wlo*t) .* filter(bcentre, 1, rx .* exp(-j*wlo*t));
  sam3 = exp(-j*wlo*t) .* filter(bcentre, 1, rx .* exp( j*wlo*t));

  % Do the math on phases using rect math

  two_phi_est_rect = conj(sam2.*sam3) .* (sam1.^2);
  phi_est = angle(two_phi_est_rect)/2;

  phi_est_av = angle(mean(two_phi_est_rect))/2;

  printf("phi_est : %3.2f rads %3.1f degrees\n", phi_est_av, phi_est_av*180/pi);

  % some plots ....

  figs = 0x1 + 0x2 + 0x20;

  if bitand(figs, 0x1)
    figure(1);
    clf;
    plot((1:N)*Fs/N, 20*log10(abs(Rx)),'markersize', 10, 'linewidth', 2);
    axis([1 Fs -80 0])
    title('Rx signal at SDR input');
  end

  if bitand(figs, 0x2)
    figure(2)
    clf
    Sam1 = (1/N)*fft(sam1,N);
    Sam2 = (1/N)*fft(sam2,N);
    Sam3 = (1/N)*fft(sam3,N);
    subplot(311)
    plot((1:N)*Fs/N, 20*log10(abs(Sam1)),'markersize', 10, 'linewidth', 2);
    axis([f-bw/2 f+bw/2 -60 0])
    title('Band Pass filtered signals');
    subplot(312)
    plot((1:N)*Fs/N, 20*log10(abs(Sam2)),'markersize', 10, 'linewidth', 2);
    axis([f+flo-bw/2 f+flo+bw/2 -60 0])
    subplot(313)
    plot((1:N)*Fs/N, 20*log10(abs(Sam3)),'markersize', 10, 'linewidth', 2);
    axis([f-flo-bw/2 f-flo+bw/2 -60 0])
  end

  if bitand(figs, 0x4)
    figure(3)
    Np = 1000;
    subplot(311)
    plot(real(sam1(1:Np)))
    title('Band Pass filtered signals');
    subplot(312)
    plot(real(sam2(1:Np)))
    subplot(313)
    plot(real(sam3(1:Np)))
  end

  if bitand(figs, 0x8)
    figure(4)
    clf;
    hist2d([real(two_phi_est_rect)' imag(two_phi_est_rect)'],50)
    xmax = 2*mean(abs(two_phi_est_rect));
    axis([-xmax xmax -xmax xmax])
    title('Scatter Plot');
  end

  if bitand(figs, 0x10)
    figure(5)
    clf;
    subplot(211)
    plot(real(two_phi_est_rect(1:Np)))
    subplot(212)
    plot(imag(two_phi_est_rect(1:Np)))
    title('two phi est rect')
  end

  if bitand(figs, 0x20)
    figure(6)
    clf
    [r theta] = hist([angle(two_phi_est_rect)/2 -pi/2 pi/2],100);
    polar([theta pi+theta],[r r])
  end

  if bitand(figs, 0x40)
    figure(7)
    clf
    plot(phi_est)
    axis([0 length(phi_est(1:Np)) -pi/2 pi/2]);
    title('phi est')
  end

  % when sampling hit any key to finish, finish after one pass in other modes

  finished = 1;
  if strcmp(mode, "sample")
    finished = length(kbhit(1));
  end
end
