% gmsk.m
% David Rowe Dec 2014
%
% GMSK modem simulation
%
% [X] plot eye diagram
% [ ] BER curves with reas match to theoretical
% [ ] spectrum - will it pass thru a HT?

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
  gmsk_states.Fs = 48000;
  gmsk_states.Rs = 4800;
  M = gmsk_states.M = gmsk_states.Fs/gmsk_states.Rs;
  global gmsk_mod_coeff;
  global gmsk_demod_coeff;

  fm_states.Fs = gmsk_states.Fs;  
  fm_states.fc = 0;  
  fm_max = fm_states.fm_max = 2400;
  fd = fm_states.fd = 1200;
  fm_states.Ts = gmsk_states.M;  
  fm_states.pre_emp = fm_states.de_emp = 0;
  fm_states.output_filter = 1;
  gmsk_states.fm_states = analog_fm_init(fm_states);

  [x i_mod] = max(gmsk_mod_coeff);
  [x i_demod] = max(gmsk_demod_coeff);
  if strcmp(gmsk_states.rx_filter,"lowpass")
    %gmsk_states.filter_delay = i_mod + i_demod;
    gmsk_states.filter_delay = i_mod + 100+21;
  elseif strcmp(gmsk_states.rx_filter,"ml")
    gmsk_states.filter_delay = i_mod + 100+35;
  else
    printf("filter type not known");
  end

  gmsk_states.Toff = 2;
  gmsk_states.dsam = dsam = gmsk_states.filter_delay;
  gmsk_states.dsym = floor(dsam/gmsk_states.M);

  % Max Likelihood 3 symbol filter
  % Matched rx filter of all possible 3 bit sequences, an attempt to
  % account for ISI.  Note filtering is in phase domain to model use
  % of legacy FM radios where FM mod.demod is not in DSP.

  ml_bits = [ 0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
  ml_symbols = zeros(8,7*M);
  for r=1:8
    for c=1:3
      ml_symbols(r,1+(c-1)*M:c*M) = -1 + 2*ml_bits(r,c);
    end
    ml_filt(r,:) = filter(gmsk_mod_coeff,1,ml_symbols(r,:));
  end
  figure(5)
  subplot(211)
  plot(ml_symbols')  
  subplot(212)
  plot(ml_filt')

  gmsk_states.ml_filt = ml_filt;
  gmsk_states.ml_bits = ml_bits;
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

  % work out delays of filter to align bits
  % plot eye diagrams, BERcurves, theoretical results, spectrum - will it pass thru a HT?
endfunction


function [rx_bits rx_out] = gmsk_demod(gmsk_states, rx)
  M = gmsk_states.M;
  nsam = length(rx);
  nsym = floor(nsam/M);
  global gmsk_demod_coeff;
  global gmsk_mod_coeff;
  wd = 2*pi*gmsk_states.fm_states.fd/gmsk_states.Fs;

  % wide filter that introduces no ISI but limits noise
  % at input to FM demod

  %fc = (4800)/(gmsk_states.Fs/2);
  %bin  = firls(200,[0 fc*(1-0.05) fc*(1+0.05) 1],[1 1 0.01 0.01]);
  %rx_filt = filter(bin, 1, rx);
  g = raised_cosine_filter(0.5,M);
  rx_filt = filter(bin, 1, rx);

  % FM demod

  rx_diff = [ 1 rx_filt(2:nsam) .* conj(rx_filt(1:nsam-1))];
  rx_out = (1/wd)*atan2(imag(rx_diff),real(rx_diff));

  % ML detector, bank of 8 filters that we run in parallel
  rx_ml_out = zeros(8,nsam);
  for r=1:8
    rx_ml_out(r,:) = filter(gmsk_states.ml_filt(r,:), 1, rx_out);
  end

  figure(6)
  clf
  nplot = 10;
  dsam = gmsk_states.dsam;
  Toff = gmsk_states.Toff;
  %subplot(8,1,1)
  st = 1+dsam+Toff;
  en = st + nplot*M-1;
  mesh(1:nplot*M, 1:8, rx_ml_out(:,st:en))
  %hold on;
  %for r=2:8
  %  %subplot(8,1,r)
  %  plot(rx_ml_out(r,1+dsam+Toff:+dsam+Toff+nplot*M))
  %end
  %hold off;

  % choose maxima

  Toff = 15;
  for i=1:nsym - (1 + dsam + Toff)/M
    s = 1 + dsam + Toff + (i-1)*M;
    [m(i) r] = max(rx_ml_out(:,s));
    rx_bits(i) = gmsk_states.ml_bits(r,2);
  end
  nplot = 200;
  figure(7)
  subplot(211)
  stem(m(1:nplot))
  subplot(212)
  stem(rx_bits(1:nplot))

  Toff = gmsk_states.Toff;
  dsam = gmsk_states.filter_delay;
  %rx_bits = rx_out(1+dsam+Toff:M:length(rx_out)) > 0;
endfunction


function sim_out = gmsk_test(sim_in)
  nsym =  sim_in.nsym;
  EbNodB = sim_in.EbNodB;
  verbose = sim_in.verbose;

  gmsk_states.rx_filter = sim_in.filter;
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

    %tx_bits = round(rand(1, nsym));
    tx_bits = ones(1, nsym);
    tx_bits(2:2:nsym) = 0;
    [tx tx_filt tx_symbols] = gmsk_mod(gmsk_states, tx_bits);
    nsam = length(tx);
    
    noise = sqrt(variance/2)*(randn(1,nsam) + j*randn(1,nsam));
    rx    = tx + noise;
    [rx_bits rx_out] = gmsk_demod(gmsk_states, rx);
      
    %tx_bits(1:10)
    %rx_bits(1:10)

    l = length(rx_bits);
    error_positions = xor(rx_bits(1:l), tx_bits(1:l));
    Nerrs = sum(error_positions);
    TERvec(ne) = Nerrs;
    BERvec(ne) = Nerrs/l;
    
    if verbose > 0
      printf("EbNo dB: %3.1f Nerrs: %d BER: %f BER Theory: %f\n", aEbNodB, Nerrs, BERvec(ne), 0.5*erfc(sqrt(0.75*EbNo)));
    end

    if verbose > 1

      figure(1)
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

      figure(2);
      nplot = 16;
      clf;
      subplot(211)
      stem(tx_symbols(1:M:nplot*M))
      axis([1 nplot -1 1])
      title('tx symbols');
      subplot(212)
      stem(rx_out(dsam+1+Toff:M:dsam+nplot*M+Toff))
      axis([1 nplot -1 1])
      title('rx symbols');

      figure(3);
      clf;
      subplot(211)
      plot(tx_filt(21+1:21+1+nplot*M))
      axis([1 nplot*M -1 1]);
      title('tx after guassian filter');
      subplot(212)
      plot(rx_out(dsam+1+Toff:dsam+nplot*M+Toff))
      axis([1 nplot*M -1 1]);
      title('rx after before sampling');

      figure(4);
      clf
      subplot(211);
      f = fft(rx);
      Tx = 20*log10(abs(f));
      plot(Tx)
      grid;
      title('GMSK Demodulator Input Spectrum');
      axis([1 10000 0 80])

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
      axis([1 10000 0 max(cs)])

      printf("Bfm: %4.0fHz %3.0f%% power bandwidth %4.0fHz = %3.2f*Rb\n", Bfm, x*100, bw, bw/Rs);
    end
  end

  sim_out.TERvec = TERvec;
  sim_out.BERvec = BERvec;
  sim_out.Rs = gmsk_states.Rs;
endfunction


function run_gmsk_single
  sim_in.filter = "lowpass";
  sim_in.nsym = 4800;
  sim_in.EbNodB = 10;
  sim_in.verbose = 2;

  sim_out = gmsk_test(sim_in);
endfunction


function run_gmsk_curves
  sim_in.nsym = 4800;
  sim_in.EbNodB = 5:10;
  sim_in.verbose = 1;

  gmsk_sim = gmsk_test(sim_in);
  Rs = gmsk_sim.Rs;
  EbNo  = 10 .^ (sim_in.EbNodB/10);
  alpha = 0.75;
  gmsk_theory.BERvec = 0.5*erfc(sqrt(alpha*EbNo));

  % BER v Eb/No curves

  figure(1); 
  clf;
  semilogy(sim_in.EbNodB, gmsk_theory.BERvec,'r;GMSK theory;')
  hold on;
  semilogy(sim_in.EbNodB, gmsk_sim.BERvec,'g;GMSK sim;')
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
  semilogy(sim_in.EbNodB+RsOnB_dB, gmsk_sim.BERvec,'g;GMSK sim;')
  hold off;
  grid("minor");
  axis([min(sim_in.EbNodB+RsOnB_dB) max(sim_in.EbNodB+RsOnB_dB) 1E-4 1])
  legend("boxoff");
  xlabel("C/No for Rs=4800 bit/s and 1 Hz noise bandwidth (dB)");
  ylabel("Bit Error Rate (BER)")

endfunction

function run_gmsk_init
  sim_in.rx_filter = "lowpass";
  gmsk_init(sim_in);
endfunction

run_gmsk_single
%run_gmsk_curves
%run_gmsk_init


