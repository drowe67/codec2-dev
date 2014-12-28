% gmsk.m
% David Rowe Dec 2014
%
% GMSK modem simulation
%
% [X] plot eye diagram
% [X] BER curves with reas match to theoretical
% [X] fine timing estimator 
% [ ] phase/freq estimator
%     + need initial acquisition and tracking
% [ ] coarse timing estimator (sync up to known test frames)
% [ ] file read/write interface
% [ ] modify for 1200 bit/s (or any) operation, ie GMSK filter coeff generation
%     or re-sampling

% Filter coeffs From:
% https://github.com/on1arf/gmsk/blob/master/gmskmodem_codec2/API/a_dspstuff.h,
% which is in turn from Jonathan G4KLX.  The demod coeffs LPF noise

global gmsk_mod_coeff = [...
 6.455906007234699e-014, 1.037067381285011e-012, 1.444835156335346e-011,...
1.745786683011439e-010, 1.829471305298363e-009, 1.662729407135958e-008,...
1.310626978701910e-007, 8.959797186410516e-007, 5.312253663302771e-006,...
2.731624380156465e-005, 1.218217140199093e-004, 4.711833994209542e-004,...
1.580581180127418e-003, 4.598383433830095e-003, 1.160259430889949e-002,...
2.539022692626253e-002, 4.818807833062393e-002, 7.931844341164322e-002,...
1.132322945270602e-001, 1.401935338024111e-001, 1.505383695578516e-001,...
1.401935338024111e-001, 1.132322945270601e-001, 7.931844341164328e-002,...
4.818807833062393e-002, 2.539022692626253e-002, 1.160259430889949e-002,...
4.598383433830090e-003, 1.580581180127420e-003, 4.711833994209542e-004,...
1.218217140199093e-004, 2.731624380156465e-005, 5.312253663302753e-006,...
8.959797186410563e-007, 1.310626978701910e-007, 1.662729407135958e-008,...
1.829471305298363e-009, 1.745786683011426e-010, 1.444835156335356e-011,...
1.037067381285011e-012, 6.455906007234699e-014];

global gmsk_demod_coeff = [...
-0.000153959924563, 0.000000000000000, 0.000167227768379, 0.000341615513437,...
0.000513334449696, 0.000667493753523, 0.000783901543032, 0.000838293462576,...
0.000805143268199, 0.000661865814384, 0.000393913058926, -0.000000000000000,...
-0.000503471198655, -0.001079755887508, -0.001671728086040, -0.002205032425392,...
-0.002594597675000, -0.002754194565297, -0.002608210441859, -0.002104352817854,...
-0.001225654870420, 0.000000000000000, 0.001494548041184, 0.003130012785731,...
0.004735238379172, 0.006109242742194, 0.007040527007323, 0.007330850462455,...
0.006821247169795, 0.005417521811131, 0.003112202160626, -0.000000000000000,...
-0.003715739376345, -0.007727358782391, -0.011638713107503, -0.014992029537478,...
-0.017304097563429, -0.018108937286588, -0.017003180218569, -0.013689829477969,...
-0.008015928769710, 0.000000000000000, 0.010154104792614, 0.022059114281395,...
0.035162729807337, 0.048781621388364, 0.062148583345584, 0.074469032280094,...
0.084982001723750, 0.093020219991183, 0.098063819576269, 0.099782731268437,...
0.098063819576269, 0.093020219991183, 0.084982001723750, 0.074469032280094,...
0.062148583345584, 0.048781621388364, 0.035162729807337, 0.022059114281395,...
0.010154104792614, 0.000000000000000, -0.008015928769710, -0.013689829477969,...
-0.017003180218569, -0.018108937286588, -0.017304097563429, -0.014992029537478,...
-0.011638713107503, -0.007727358782391, -0.003715739376345, -0.000000000000000,...
0.003112202160626, 0.005417521811131, 0.006821247169795, 0.007330850462455,...
0.007040527007323, 0.006109242742194, 0.004735238379172, 0.003130012785731,...
0.001494548041184, 0.000000000000000, -0.001225654870420, -0.002104352817854,...
-0.002608210441859, -0.002754194565297, -0.002594597675000, -0.002205032425392,...
-0.001671728086040, -0.001079755887508, -0.000503471198655, -0.000000000000000,...
0.000393913058926, 0.000661865814384, 0.000805143268199, 0.000838293462576,...
0.000783901543032, 0.000667493753523, 0.000513334449696, 0.000341615513437,...
0.000167227768379, 0.000000000000000, -0.000153959924563];

rand('state',1); 
randn('state',1);
graphics_toolkit ("gnuplot");
fm;

function gmsk_states = gmsk_init(gmsk_states)

  % general 

  gmsk_states.Fs = 48000;
  gmsk_states.Rs = 4800;
  M = gmsk_states.M = gmsk_states.Fs/gmsk_states.Rs;
  global gmsk_mod_coeff;
  global gmsk_demod_coeff;

  % set up FM modulator

  fm_states.Fs = gmsk_states.Fs;  
  fm_states.fc = 0;  
  fm_max = fm_states.fm_max = 2400;
  fd = fm_states.fd = 1200;
  fm_states.Ts = gmsk_states.M;  
  fm_states.pre_emp = fm_states.de_emp = 0;
  fm_states.output_filter = 1;
  gmsk_states.fm_states = analog_fm_init(fm_states);

  % configure delays for demod type

  [x i_mod] = max(gmsk_mod_coeff);
  [x i_demod] = max(gmsk_demod_coeff);

  if gmsk_states.coherent_demod
    gmsk_states.filter_delay = i_mod + i_mod;
    gmsk_states.Toff = -4;
  else
    gmsk_states.filter_delay = i_mod + i_demod + 100;
    gmsk_states.Toff = 3;
 end

  gmsk_states.dsam = dsam = gmsk_states.filter_delay;

endfunction


function [tx tx_filt tx_symbols] = gmsk_mod(gmsk_states, tx_bits)
  M = gmsk_states.M;
  nsym = length(tx_bits);
  nsam = nsym*M;
  global gmsk_mod_coeff;

  % NRZ sequence of symbols

  tx_symbols = zeros(1,nsam);
  for i=1:nsym
    tx_symbols(1+(i-1)*M:i*M) = -1 + 2*tx_bits(i);
  end

  tx_filt = filter(gmsk_mod_coeff, 1, tx_symbols);
  tx = analog_fm_mod(gmsk_states.fm_states, tx_filt);
endfunction


function [rx_bits rx_int rx_filt] = gmsk_demod(gmsk_states, rx)
  M = gmsk_states.M;
  Rs = gmsk_states.Rs;
  Fs = gmsk_states.Fs;
  nsam = length(rx);
  nsym = floor(nsam/M);
  global gmsk_demod_coeff;
  global gmsk_mod_coeff;
  wd = 2*pi*gmsk_states.fm_states.fd/gmsk_states.Fs;
  Toff = gmsk_states.Toff;
  dsam = gmsk_states.filter_delay;

  if gmsk_states.coherent_demod

    % See IEEE Trans on Comms, Muroyta et al, 1981, "GSM Modulation
    % for Digital Radio Telephony" Fig 8:

    % matched filter

    rx_filt = filter(gmsk_mod_coeff, 1, rx);

    % Property of MSK that re and im arms are sequences of 2T
    % long symbols, can be demodulated like QPSK with matched filter
    % and integrate and dump.

    % integrate energy in symbols 2T long in re and im arms
    % note this could be combined with matched filter

    rx_int = conv(rx_filt,ones(1,2*M));

    % TODO implement and test phase/fine freq tracking loop

    dsam = 41 + M;

    % timing estimate todo: clean up all these magic numbers

    x = abs(real(rx_int));
    w = 2*pi*(Rs/2)/Fs;
    X = x(1:nsam) * exp(j*w*(0:nsam-1))';
    timing_angle = angle(X);
    timing_adj = timing_angle*2*M/(2*pi);

    timing_adj = 1.4; % lock for now while we play with phase

    printf("timing angle: %f adjustment: %f\n", angle(X), timing_adj);
    dsam -= floor(timing_adj) - 2;

    if gmsk_states.phase_track

      % DCO design from "Introduction To Phase-Lock Loop System Modeling", Wen Li
      % http://www.ece.ualberta.ca/~ee401/parts/data/PLLIntro.pdf

      eta = 0.707;
      wn = 2*pi*20;
      Ts = 1/Fs;
      g1 = 1 - exp(-2*eta*wn*Ts);
      g2 = 1 + exp(-2*eta*wn*Ts) - 2*exp(-eta*wn*Ts)*cos(wn*Ts*sqrt(1-eta*eta));
      Gpd = 2/pi;
      Gvco = 1;
      G1 = g1/(Gpd*Gvco);  G2 = g2/(Gpd*Gvco);
      printf("g1: %e g2: %e G1: %e G2: %e\n", g1, g2, G1, G2);

      filt_prev = dco = lower = ph_err_filt = ph_err = 0;
      dco_log = filt_log = zeros(1,nsam);

      % we need a sine wave at the timing clock frequency

      timing_clock_phase = timing_angle;
 
      for i=1:nsam
        timing_clock_phase += (2*pi)/(2*M);
        rx_int(i) *= exp(-j*dco);
        ph_err = sign(real(rx_int(i))*imag(rx_int(i)))*sin(timing_clock_phase);
        lower = ph_err*G2 + lower;
        filt  = ph_err*G1 + lower;
        dco   = dco + filt_prev;
        filt_prev = filt;
        filt_log(i) = filt;
        dco_log(i) = dco;
      end

      figure(6);
      clf
      subplot(211);
      plot(filt_log);
      subplot(212);
      plot(dco_log/pi);
      %axis([1 nsam -0.5 0.5])
    end

    % sample symbols at end of integration

    re_syms = real(rx_int(1+dsam+Toff:2*M:nsam));
    im_syms = imag(rx_int(1+dsam+M+Toff:2*M:nsam));

    % XORs/adders on the RHS of Muroyta et al Fig 8 (a) and (b).  We
    % simulate digital logic bit stream at clock rate Rs, even though
    % we sample integrators at rate Rs/2.  I can't explain how and why
    % this logic works/is required.  I think it can be worked out from
    % comparing MSK demods.

    l = length(re_syms);
    l2 = 2*l;
    re_bits = zeros(1,l2);
    im_bits = zeros(1,l2);
    clk_bits = zeros(1,l2);
    for i=1:l-1
      re_bits(2*(i-1)+1)  = re_syms(i) > 0;
      re_bits(2*(i-1)+2)  = re_syms(i) > 0;
      im_bits(2*(i-1)+2)  = im_syms(i) > 0;
      im_bits(2*(i-1)+3)  = im_syms(i) > 0;
      clk_bits(2*(i-1)+1) = 0;
      clk_bits(2*(i-1)+2) = 1;
    end

    rx_bits = bitxor(bitxor(re_bits,im_bits),  clk_bits);
    rx_bits = rx_bits(2:length(rx_bits)-1);
  else
    % non-coherent demod

    % filter to get rid of most of noise before FM demod, but doesnt
    % introduce any ISI

    fc = (4800)/(gmsk_states.Fs/2);
    bin  = firls(200,[0 fc*(1-0.05) fc*(1+0.05) 1],[1 1 0.01 0.01]);
    rx_filt = filter(bin, 1, rx);

    % FM demod

    rx_diff = [ 1 rx_filt(2:nsam) .* conj(rx_filt(1:nsam-1))];
    rx_out = (1/wd)*atan2(imag(rx_diff),real(rx_diff));

    % low pass filter, trade off betwen ISI and removing noise

    rx_out = filter(gmsk_demod_coeff, 1, rx_out);
    
    rx_bits = real(rx_out(1+dsam+Toff:M:length(rx_out)) > 0);
  end

endfunction


function [freq_offset_est phase_offset_est] = gmsk_est_freq_phase(gmsk_states, rx, verbose)
  Fs = gmsk_states.Fs;
  Rs = gmsk_states.Rs;

  % one frame should have preamble at start and end, enough to see in DFT 

  f = fft(rx);
  nsam = length(rx);

  % Look for line at Rs/2.  Just look on +ve side for now, a more
  % reliable method would be to check for mirror image of line at -ve,
  % look for line in range +/- 100 Hz
  
  start_bin = floor((Rs/2 - 100)*nsam/Fs) + 1;
  stop_bin = floor((Rs/2 + 100)*nsam/Fs) + 1;
  [max_val max_bin] = max(abs(f(start_bin:stop_bin)));
  printf("nsam: %d start_bin: %d start_bin: %d max_bin: %d\n", nsam, start_bin, stop_bin, max_bin);

  % map max_bin to frequency offset

  freq_offset_est = (max_bin - floor(100*nsam/Fs) - 1)*Fs/nsam;  
  phase_offset_est = angle(f(max_bin)) - 0.935153; % magic number due to initial filter delay?

  if verbose
    printf("freq_offset_est: %f phase_offset_est: %f\n", freq_offset_est, phase_offset_est);
    figure(6)
    plot((1:nsam)*Fs/nsam, 20*log10(abs(f)));
    axis([1 Rs 0 80]);
  end

endfunction


function sim_out = gmsk_test(sim_in)
  nsym =  sim_in.nsym;
  EbNodB = sim_in.EbNodB;
  verbose = sim_in.verbose;

  gmsk_states.coherent_demod = sim_in.coherent_demod;
  gmsk_states.phase_track = 0;
  gmsk_states = gmsk_init(gmsk_states);
  M = gmsk_states.M;
  Fs = gmsk_states.Fs;
  Rs = gmsk_states.Rs;
  Bfm = gmsk_states.fm_states.Bfm;
  dsam = gmsk_states.dsam;
  Toff = gmsk_states.Toff;
 
  for ne = 1:length(EbNodB)
    aEbNodB = EbNodB(ne);
    EbNo = 10^(aEbNodB/10);
    variance = Fs/(Rs*EbNo);

    tx_bits = round(rand(1, nsym));
    tx_bits = ones(1, nsym);
    %tx_bits = zeros(1, nsym);
    tx_bits(1:2:nsym) = 0;
    [tx tx_filt tx_symbols] = gmsk_mod(gmsk_states, tx_bits);
    nsam = length(tx);
    %tx_bits(1:10)
    noise = sqrt(variance/2)*(randn(1,nsam) + j*randn(1,nsam));
    rx    = tx*exp(j*pi/2) + noise;

    [rx_bits rx_out rx_filt] = gmsk_demod(gmsk_states, rx(1:length(rx)));
      
    l = length(rx_bits);
    error_positions = xor(rx_bits(1:l), tx_bits(1:l));
    Nerrs = sum(error_positions);
    TERvec(ne) = Nerrs;
    BERvec(ne) = Nerrs/l;
    
    if verbose > 0
      printf("EbNo dB: %3.1f Nerrs: %d BER: %f BER Theory: %f\n", aEbNodB, Nerrs, BERvec(ne), 0.5*erfc(sqrt(0.75*EbNo)));
    end

    if verbose > 1

      if gmsk_states.coherent_demod == 0
        figure(1)
        clf
        eyesyms = 2;
        plot(rx_out(dsam+1+Toff:dsam+eyesyms*M+Toff))
        hold on;
        for i=1:10
          st = dsam+1+Toff+i*eyesyms*M;
          en = st + eyesyms*M;
          plot(rx_out(st:en))
        end
        hold off;
        axis([0 eyesyms*M -2 2]);
        title('Eye Diagram');
      else
        figure(1);
        nplot = 16;
        clf;
        subplot(211)
        plot(real(rx_filt(1:nplot*M)))
        axis([1 nplot*M -1 1])
        title('Matched Filter');
        subplot(212)
        plot(imag(rx_filt(1:nplot*M)))
        axis([1 nplot*M -1 1])

        figure(2);
        nplot = 16;
        clf;
        subplot(211)
        plot(real(rx_out(1:nplot*M))/(2*M))
        title('Integrator');
        axis([1 nplot*M -1 1])
        subplot(212)
        plot(imag(rx_out(1:nplot*M)/(2*M)))
        axis([1 nplot*M -1 1])
     end

      figure(3)
      clf
      subplot(211)
      stem(tx_bits(1:20))
      title('Tx Bits')
      subplot(212)
      stem(rx_bits(1:20))
      title('Rx Bits')

      figure(4);
      clf
      subplot(211);
      f = fft(rx);
      Tx = 20*log10(abs(f));
      plot(Tx)
      grid;
      title('GMSK Demodulator Input Spectrum');
      axis([1 5000 0 80])

      subplot(212)
      f = fft(tx);
      f = f(1:length(f)/2);
      cs = cumsum(abs(f).^2);
      plot(cs)
      hold on;
      x = 0.99;
      tots = x*sum(abs(f).^2);
      xpercent_pwr = find(cs > tots);
      bw = 2*xpercent_pwr(1);
      plot([1 Fs/2],[tots tots],'r')
      plot([bw/2 bw/2],[0 tots],'r')
      hold off;  
      title("Cumulative Power");
      grid;
      axis([1 5000 0 max(cs)])

      printf("Bfm: %4.0fHz %3.0f%% power bandwidth %4.0fHz = %3.2f*Rb\n", Bfm, x*100, bw, bw/Rs);

    end
  end

  sim_out.TERvec = TERvec;
  sim_out.BERvec = BERvec;
  sim_out.Rs = gmsk_states.Rs;
endfunction


function run_gmsk_single
  sim_in.coherent_demod = 1;
  sim_in.nsym = 480;
  sim_in.EbNodB = 6;
  sim_in.verbose = 2;

  sim_out = gmsk_test(sim_in);
endfunction


function run_gmsk_curves
  sim_in.coherent_demod = 1;
  sim_in.nsym = 4800;
  sim_in.EbNodB = 2:10;
  sim_in.verbose = 1;

  gmsk_coh = gmsk_test(sim_in);

  sim_in.coherent_demod = 0;
  gmsk_noncoh = gmsk_test(sim_in);

  Rs = gmsk_coh.Rs;
  EbNo  = 10 .^ (sim_in.EbNodB/10);
  alpha = 0.75; % guess for BT=0.5 GMSK
  gmsk_theory.BERvec = 0.5*erfc(sqrt(alpha*EbNo));

  % BER v Eb/No curves

  figure(1); 
  clf;
  semilogy(sim_in.EbNodB, gmsk_theory.BERvec,'r;GMSK theory;')
  hold on;
  semilogy(sim_in.EbNodB, gmsk_coh.BERvec,'g;GMSK sim coherent;')
  semilogy(sim_in.EbNodB, gmsk_noncoh.BERvec,'b;GMSK sim non-coherent;')
  hold off;
  grid("minor");
  axis([min(sim_in.EbNodB) max(sim_in.EbNodB) 1E-4 1])
  legend("boxoff");
  xlabel("Eb/No (dB)");
  ylabel("Bit Error Rate (BER)")

  % BER v C/No (1 Hz noise BW and Eb=C/Rs=1/Rs)
  % Eb/No = (C/Rs)/(1/(N/B))
  % C/N   = (Eb/No)*(Rs/B)

  RsOnB_dB = 10*log10(Rs/1);
  figure(2); 
  clf;
  semilogy(sim_in.EbNodB+RsOnB_dB, gmsk_theory.BERvec,'r;GMSK theory;')
  hold on;
  semilogy(sim_in.EbNodB+RsOnB_dB, gmsk_coh.BERvec,'g;GMSK sim coherent;')
  semilogy(sim_in.EbNodB+RsOnB_dB, gmsk_noncoh.BERvec,'b;GMSK sim non-coherent;')
  hold off;
  grid("minor");
  axis([min(sim_in.EbNodB+RsOnB_dB) max(sim_in.EbNodB+RsOnB_dB) 1E-4 1])
  legend("boxoff");
  xlabel("C/No for Rs=4800 bit/s and 1 Hz noise bandwidth (dB)");
  ylabel("Bit Error Rate (BER)")

endfunction

% TODO: [ ] test over range of freq, phase, coarse timing and Eb/No values
%       [ ] Modify for Rs=1200, e.g. mod filter above
%       [ ] extra stuff at begining/end for filter delays
   
function run_test_freq_offset
  verbose = 1;
  aEbNodB = 6;
  phase_offset = 0;
  freq_offset  = 1;
  nsym = 4800;

  gmsk_states.coherent_demod = 1;
  gmsk_states.phase_track    = 1;
  gmsk_states = gmsk_init(gmsk_states);
  Fs = gmsk_states.Fs;
  Rs = gmsk_states.Rs;

  % A frame consists of nsym random data bits.  Some experimentation
  % has shown they need to be random for the sync algorithms to work

  tx_bits = round(rand(1, nsym));

  [tx tx_filt tx_symbols] = gmsk_mod(gmsk_states, tx_bits);
  nsam = length(tx);

  EbNo = 10^(aEbNodB/10);
  variance = Fs/(Rs*EbNo);
  noise = sqrt(variance/2)*(randn(1,nsam) + j*randn(1,nsam));
  w  = (0:nsam-1)*2*pi*freq_offset/Fs + phase_offset;
  figure(2)
  plot(w/pi)

  rx = tx.*exp(j*w) + noise;

  [freq_offset_est phase_offset_est] = gmsk_est_freq_phase(gmsk_states, rx, verbose);
  w_est  = (0:nsam-1)*2*pi*freq_offset_est/Fs + phase_offset_est;
  %rx = rx.*exp(-j*w_est);

  [rx_bits rx_out rx_filt] = gmsk_demod(gmsk_states, rx);
  l = length(rx_bits);
  error_positions = xor(rx_bits(1:l), tx_bits(1:l));
  Nerrs = sum(error_positions);
  ber = Nerrs/l;
  printf("Eb/No: %3.1f freq_offset: %4.1f phase_offset: %4.3f Nerrs: %d BER: %f\n", 
         aEbNodB, freq_offset, phase_offset, Nerrs, ber);
  figure(1)
  clf
  subplot(311)
  stem(real(tx_bits(500:520)))
  subplot(312)
  stem(real(rx_bits(500:520)))
  subplot(313)
  stem(real(cumsum(error_positions)))
endfunction
    
%run_gmsk_single
%run_gmsk_curves
%run_gmsk_init
run_test_freq_offset

