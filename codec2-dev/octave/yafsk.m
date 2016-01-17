% yafsk.m
% Yet-another-FSK
% Brady O'Brien 20 October 2015
%
% This is a model for the first attempt at a C FSK modem. Based on fsk_horus and maybe a little fsk4.
% First revision will just be 2fsk
% [x] - Modem init and struct def
% [x] - Direct SDR modulator, probably not FM based
% [x] - Direct SDR non-coherent demodulator
%    [x] - Core demodulation routine
%    |o| - Timing offset estimation
%    < >- Freq. offset estimation
%    (+) - Bit slip, maybe
% { } - Port sim from fsk_horus
% [ ] - The C port
% [ ] - Some stuff to verify the C port
% [ ] - 4FSK variant
%    ( ) - All of that other stuff, but for 4fsk
%clear all;

graphics_toolkit('gnuplot');
%fm

pkg load signal;

%Basic parameters for a simple FSK modem
fsk_setup_info.Rs = 4800;  % Symbol rate
fsk_setup_info.nfsk = 2;      % Number of unique symbols. Must be 2.
fsk_setup_info.P = 5;		%Something something fine timing est
fsk_setup_info.Fs = 48000; % Sample frequency
fsk_setup_info.timing_syms = 10; %How many symbols over which to figure fine timing
fsk_setup_info.Fsym = fsk_setup_info.Rs; %Symbol spacing
fsk_setup_info.txmap = @(bits) bits+1; %Map TX bits to 2fsk symbols
fsk_setup_info.rxmap = @(syms) syms==1; %Map 2fsk RX symbols to bits
global fsk_setup = fsk_setup_info

function states = yafsk_init(fsk_config)
  Fs = states.Fs = fsk_config.Fs;
  Rs = states.Rs = fsk_config.Rs;
  nfsk = states.nfsk = fsk_config.nfsk;
  Ts = states.Ts = Fs/Rs;
  Fsym = states.Fsym = fsk_config.Fsym;
  states.config = fsk_config;
  P = states.P = fsk_config.P;
  timing_syms = states.timing_syms = fsk_config.timing_syms;
  
  if nfsk != 2
    error("Gotta be 2fsk")
  endif
 
  %Symbol frequencies. Fixed to intervals of Fsym.
  states.fsyms = [-(Fsym/2) (Fsym/2)];
  states.tx_phase = 0;

  states.dc = zeros(1,nfsk);
  states.rx_phi = ones(1,nfsk);
  states.isamp = 0;
  states.ssamp = 0;
  states.sums = zeros(1,nfsk);
  states.ssums = zeros(1,nfsk);
  timing_db1 = timing_db2 = zeros(1,timing_syms*(Ts/P))
  states.timing_db1 = timing_db1;
  states.timing_db2 = timing_db2;
endfunction

function [tx states] = yafsk_mod(states,bits)
  Ts = states.Ts;
  Fs = states.Fs;
  fsyms = states.fsyms;
  tx_phase = states.tx_phase;
  %Map bits into symbols
  syms = states.config.txmap(bits);
  tx = zeros(1,Ts*length(syms));
  
  for ii = (1:length(syms))
    cur_sym_f = fsyms(syms(ii));
    tx_phase_i = tx_phase;
    for jj = (1:Ts)
        tx_phase_i = tx_phase + jj*2*pi*cur_sym_f/Fs;
        tx((ii-1)*Ts+jj) = exp(j*tx_phase_i);   
    end
    tx_phase = tx_phase + Ts*2*pi*cur_sym_f/Fs;
    if tx_phase>2*pi
        tx_phase = tx_phase-2*pi;
    elseif tx_phase<-2*pi
        tx_phase = tx_phase+2*pi;
    endif
    %tx_phase_vec = tx_phase + (1:Ts)*2*pi*cur_sym_f/Fs;
    %tx((ii-1)*Ts+1:ii*Ts) = exp(j*tx_phase_vec);
  end
  states.tx_phase = tx_phase;
endfunction

function d = idmp(data, M)
    d = zeros(1,length(data)/M);
    for i = 1:length(d)
      d(i) = sum(data(1+(i-1)*M:i*M));
    end
endfunction

function [bits states phis softsyms] = yafsk_demod_2a(states,rx)
  fine_timing = 1;
  Fs = states.Fs;
  Rs = states.Rs;
  Ts = states.Ts;
  nfsk = states.nfsk;
  P = 5;


  phy_f1 = states.rx_phi(1);
  phy_f2 = states.rx_phi(2);

  dphase_f1 = exp(states.fsyms(1)*-j*2*pi/Fs);
  dphase_f2 = exp(states.fsyms(2)*-j*2*pi/Fs);

  sum_f1 = states.sums(1);
  sum_f2 = states.sums(2);
  
  ssum_f1 = states.ssums(1);
  ssum_f2 = states.ssums(2);
  
  timing_db1 = states.timing_db1;
  timing_db2 = states.timing_db2;
  
  ssamp = states.ssamp;
  isamp = states.isamp;
  
  symcnt = 1;
  subcnt = 1;
  syms = [0];
  softsyms = [0];
  sums1 = [0];
  sums2 = [0];
  phis = [0];
  isamp = 1;
  isub = 1;

  ssum_f1 = 0;
  ssum_f2 = 0;
  ssamp=0;
 
  timing_syms = states.timing_syms;
  timing_nudge = .09; %How far to 'nudge' the sampling point
                      %This really ought to be fixed somewhere else
                      
  timing_samps = timing_syms*(Ts/P);

  for ii = (1:length(rx))
    phy_f1 *= dphase_f1;   %Spin the oscillators
    phy_f2 *= dphase_f2;

    dcs_f1 = rx(ii)*phy_f1; %Figure out the DC
    dcs_f2 = rx(ii)*phy_f2;

    sum_f1 += dcs_f1; %Integrate
    sum_f2 += dcs_f2;
    ssum_f1 += dcs_f1;
    ssum_f2 += dcs_f2;


    %Frequency of timing tracking nonlinearity
    w = (2*pi*Rs)/(Rs*P);

    %increment symbol and timing sub-loop counters
    ssamp += 1;
    isamp += 1;
    if isamp>=Ts %If it's time to take a sample and spit out a symbol..
      syms(symcnt) = (abs(sum_f1)>abs(sum_f2))+1; %Spit out a symbol
      softsyms(symcnt) = abs(sum_f1) - abs(sum_f2);
      symcnt += 1;
      
      %Fine timing estimation and adjustment
      %Re-integrate over entire symbol peiod and apply nonlinearity
      %To-do - rewrite in more c-able loop
      timing_phi = 0;
      for jj=(1:(length(timing_db1)-P))
	    f1 = sum(timing_db1(jj:jj+P)); 
        f2 = sum(timing_db2(jj:jj+P));
        tdmd = (abs(f1)-abs(f2))^2; 
        timing_phi += tdmd*exp(-j*w*jj);
      end
      
      %Take angle of nonlinear spectral line
      timing_phi = angle(timing_phi);
      %get another number
      norm_phi = timing_phi/(2*pi);
      norm_timing = norm_phi*P;
      
      %Ideal fine timing point, determined experimentally
      norm_tgt = -1.37;
      
      %move sampling point a bit forward or backward
      if(norm_timing>norm_tgt && norm_timing<norm_tgt+2.5)
          isamp += timing_nudge;
          ssamp += timing_nudge;
      else
          isamp -= timing_nudge;
          ssamp -= timing_nudge;
      endif
      
      %save timing for debugging
      phis(symcnt) = norm_timing;
      
      sum_f1 = 0;    %Reset integrators
      sum_f2 = 0;
      isamp -= Ts;    %Reset integrator count
     
      if(mod(symcnt,10000)==0)
        ab_f1 = abs(phy_f1)
        phy_f1 = phy_f1/ab_f1;
        ab_f2 = abs(phy_f2)
        phy_f2 = phy_f2/ab_f2;
	  endif

    endif

    if ssamp>= (Ts/P)
    
      %save timing samples
      timing_db1(1:timing_samps-1)=timing_db1(2:timing_samps);
      timing_db1(timing_samps) = ssum_f1;
      timing_db2(1:timing_samps-1)=timing_db2(2:timing_samps);
      timing_db2(timing_samps) = ssum_f2;


      %Reset integrators and sampling counter
      ssum_f1 = 0; 
      ssum_f2 = 0;
      ssamp -= Ts/P; 
    endif

  end

  states.rx_phy(1) = phy_f1;
  states.rx_phy(2) = phy_f2;
 
  states.sums(1) = sum_f1;
  states.sums(2) = sum_f2;

  states.ssum(1) = ssum_f1;
  states.ssum(2) = ssum_f2;

  states.ssamp = ssamp;
  states.isamp = isamp;
  
  states.timing_db1 = timing_db1;
  states.timing_db2 = timing_db2;

  bits = states.config.rxmap(syms);
  
  
endfunction


function [bits states] = yafsk_demod_2(states,rx)
  fine_timing = 1;
  Fs = states.Fs;
  Rs = states.Rs;
  Ts = states.Ts;
  nfsk = states.nfsk;

  rx = rx(fine_timing:length(rx));
  sym_phases = (1:length(rx)).*rot90(states.fsyms)*2*pi/Fs;

  sym_mixers = exp(-j*sym_phases);
  rx_mixd = repmat(rx,nfsk,1).*sym_mixers;
  
  dc1 = abs(idmp(rx_mixd(1,1:length(rx_mixd)),Ts));
  dc2 = abs(idmp(rx_mixd(2,1:length(rx_mixd)),Ts));

  t=(1:length(dc1));
  
  plot(t,dc1,t,dc2)

endfunction

% Simulation thing, shamelessly taken from fsk_horus.m
% simulation of tx and rx side, add noise, channel impairments ----------------------

function run_sim
  global fsk_setup;
  frames = 10;
  EbNodB = 26;
  timing_offset = 0.0; % see resample() for clock offset below
  test_frame_mode = 2;
  fading = 0;          % modulates tx power at 2Hz with 20dB fade depth, 
                       % to simulate balloon rotating at end of mission
  df     = 0;          % tx tone freq drift in Hz/s

  more off
  rand('state',100); 
  randn('state',100);
  states = yafsk_init(fsk_setup);
  N = states.N;
  %P = states.P;
  Rs = states.Rs;
  nsym = states.nsym = 4096;
  Fs = states.Fs;
  states.df = df;

  EbNo = 10^(EbNodB/10);
  variance = states.Fs/(states.Rs*EbNo);

  % set up tx signal with payload bits based on test mode

  if test_frame_mode == 1
     % test frame of bits, which we repeat for convenience when BER testing
    test_frame = round(rand(1, states.nsym));
    tx_bits = [];
    for i=1:frames+1
      tx_bits = [tx_bits test_frame];
    end
  end
  if test_frame_mode == 2
    % random bits, just to make sure sync algs work on random data
    tx_bits = round(rand(1, states.nsym*(frames+1)));
  end
  if test_frame_mode == 3
    % ...10101... sequence
    tx_bits = zeros(1, states.nsym*(frames+1));
    tx_bits(1:2:length(tx_bits)) = 1;
  end

  if test_frame_mode == 4

    % load up a horus msg from disk and modulate that

    test_frame = load("horus_msg.txt");
    ltf = length(test_frame);
    ntest_frames = ceil((frames+1)*nsym/ltf);
    tx_bits = [];
    for i=1:ntest_frames
      tx_bits = [tx_bits test_frame];
    end
  end

  tx = yafsk_mod(states, tx_bits);

  %tx = resample(tx, 1000, 1001); % simulated 1000ppm sample clock offset

  if fading
     ltx = length(tx);
     tx = tx .* (1.1 + cos(2*pi*2*(0:ltx-1)/Fs))'; % min amplitude 0.1, -20dB fade, max 3dB
  end

  noise = sqrt(variance)*randn(length(tx),1);
  rx    = tx + noise;
  
  % dump simulated rx file
  ftx=fopen("ya_fsk_rx.raw","wb"); rxg = rx*1000; fwrite(ftx, rxg, "short"); fclose(ftx);

  timing_offset_samples = round(timing_offset*states.Ts);
  st = 1 + timing_offset_samples;
  rx_bits_buf = zeros(1,2*nsym);
  Terrs = Tbits = 0;
  state = 0;
  x_log = [];
  norm_rx_timing_log = [];
  nerr_log = [];
  f1_int_resample_log = [];
  f2_int_resample_log = [];
  f1_log = f2_log = [];
  EbNodB_log = [];
  rx_bits_log = [];
  rx_bits_sd_log = [];

  for f=1:frames

    % extract nin samples from input stream

    nin = states.nin;
    en = st + states.nin - 1;
    sf = rx(st:en);
    st += nin;

    % demodulate to stream of bits

    [rx_bits states] = yafsk_demod(states, sf);
    rx_bits_buf(1:nsym) = rx_bits_buf(nsym+1:2*nsym);
    rx_bits_buf(nsym+1:2*nsym) = rx_bits;
    rx_bits_log = [rx_bits_log rx_bits];
    rx_bits_sd_log = [rx_bits_sd_log states.rx_bits_sd];

    norm_rx_timing_log = [norm_rx_timing_log states.norm_rx_timing];
    x_log = [x_log states.x];
    f1_int_resample_log = [f1_int_resample_log abs(states.f1_int_resample)];
    f2_int_resample_log = [f2_int_resample_log abs(states.f2_int_resample)];
    f1_log = [f1_log states.f1];
    f2_log = [f2_log states.f2];
    EbNodB_log = [EbNodB_log states.EbNodB];

    % frame sync based on min BER

    if test_frame_mode == 1
      nerrs_min = nsym;
      next_state = state;
      if state == 0
        for i=1:nsym
          error_positions = xor(rx_bits_buf(i:nsym+i-1), test_frame);
          nerrs = sum(error_positions);
          if nerrs < nerrs_min
            nerrs_min = nerrs;
            coarse_offset = i;
          end
        end
        if nerrs_min < 3
          next_state = 1;
          %printf("%d %d\n", coarse_offset, nerrs_min);
        end
      end

      if state == 1  
        error_positions = xor(rx_bits_buf(coarse_offset:coarse_offset+nsym-1), test_frame);
        nerrs = sum(error_positions);
        Terrs += nerrs;
        Tbits += nsym;
        nerr_log = [nerr_log nerrs];
      end

      state = next_state;
    end
  end

  if test_frame_mode == 1
    printf("frames: %d Tbits: %d Terrs: %d BER %4.3f\n", frames, Tbits, Terrs, Terrs/Tbits);
  end

  if test_frame_mode == 4
    extract_and_print_packets(states, rx_bits_log, rx_bits_sd_log)
  end

  figure(1);
  plot(f1_int_resample_log,'+')
  hold on;
  plot(f2_int_resample_log,'g+')
  hold off;

  figure(2)
  clf
  m = max(abs(x_log));
  plot(x_log,'+')
  axis([-m m -m m])
  title('fine timing metric')

  figure(3)
  clf
  subplot(211)
  plot(norm_rx_timing_log);
  axis([1 frames -1 1])
  title('norm fine timing')
  subplot(212)
  plot(nerr_log)
  title('num bit errors each frame')

  figure(4)
  clf
  plot(real(rx(1:Fs)))
  title('rx signal at demod input')

  figure(5)
  clf
  plot(f1_log)
  hold on;
  plot(f2_log,'g');
  hold off;
  title('tone frequencies')
  axis([1 frames 0 Fs/2])

  figure(6)
  clf
  plot(EbNodB_log);
  title('Eb/No estimate')

endfunction

% Bit error rate test ----------------------------------------------------------
% Params - aEsNodB - EbNo in decibels
%        - timing_offset - how far the fine timing is offset
%        - bitcnt - how many bits to check
%        - demod_fx - demodulator function
% Returns - ber - teh measured BER
%         - thrcoh - theory BER of a coherent demod
%         - thrncoh - theory BER of non-coherent demod
function [ber thrcoh thrncoh rxphis] = nfbert_2(aEsNodB,modem_config, bitcnt=12000, timing_offset = 1, freq_offset = 0, burst = 0,samp_clk_offset = 0)

  rand('state',1); 
  randn('state',1);
  
  %How many bits should this test run?
  %bitcnt = 12000;
  
  framesync = []
  framehdr = [1 0 1 0 1 0 1 0 0 1 1 1 0 1 0 0]
  convhdr = framehdr;%[.5 -.5 .5 -.5 1 -1 1 -1 -1 1 1 1 -1 1 -1 -1];
  %framehdr = [1 0 1 0 0 0 1 0 1 0 0 1]
  test_bits = [framesync framehdr (rand(1,bitcnt)>.5)]; %Random bits. Pad with zeros to prime the filters
  states.M = 1;
  states = yafsk_init(modem_config);
  
  %Set this to 0 to cut down on the plotting
  states.verbose = 1; 
  Fs = states.Fs;
  Rb = states.Rs;  % Multiply symbol rate by 2, since we have 2 bits per symbol
  
  tx = yafsk_mod(states,test_bits);

  tx = resample(tx, 1000, 1000 + samp_clk_offset);
  
  %simulate a single frame in a pool of noise
  if(burst)
    tx = [zeros(1,Fs/2) tx zeros(1,Fs/2)];
  end


  %add noise here
  %shamelessly copied from gmsk.m
  EsNo = 10^(aEsNodB/10);
  EbNo = EsNo
  variance = Fs/(Rb*EbNo);
  nsam = length(tx);
  noise = sqrt(variance/2)*(randn(1,nsam) + j*randn(1,nsam));
  rx    = tx*exp(j*pi/2) + noise;
  t = (1:length(rx));
  rx    = rx .* exp(t*2*pi*j*(freq_offset/Fs));
  
  rx    = rx(timing_offset:length(rx));
  
  [rx_bits states rxphis sbits] = yafsk_demod_2a(states,rx);
  
  
  ber = 1;
  
  hsig = -1*fliplr((framehdr*2)-1);
  figure(3);
  corrfd = conv(hsig,sbits);
  for ii = (1:length(corrfd)-length(hsig))
    secte = sum(abs(sbits(ii:ii+length(hsig))).^2);
	corrfd = corrfd(ii)/secte;
  end
  %thing to account for offset from input data to output data
  %No preamble detection yet
  figure(4);
  plot(rxphis);
  ox = 1;
  
  if(burst)
    orange = (1:length(rx_bits)-length(test_bits)-1);
  else
    orange = (1:100)
  end
  for offset = orange
    if(burst)
      perr = xor(test_bits,rx_bits(offset:offset-1+length(test_bits)));
      nerr = sum(perr);
      bern = nerr / length(test_bits);
	else
	  perr = xor(test_bits(offset:length(rx_bits)),rx_bits(1:length(rx_bits)+1-offset));
      nerr = sum(perr);      
      bern = nerr/(bitcnt-offset);
    end
    
    if(bern < ber)
      ox = offset;
      best_nerr = nerr;
      xerr = perr;
    end
    ber = min([ber bern]);
  end
  
  %Try to find frame header
  for offset = (1:length(rx_bits)-length(framehdr))
	hd = sum(xor(framehdr,rx_bits(offset:offset-1+length(framehdr))));
	if(hd<=2)
		printf("Found possible header at offset %d\n",offset);
	endif
  end
  
  figure(5);
  stairs(1.0*xerr);
  
  offset = ox;
  printf("\ncoarse timing: %d nerr: %d\n", offset, best_nerr);

  % Coherent BER theory
  thrcoh = .5*erfc(sqrt(.5*EbNo));

  % non-coherent BER theory calculation
  % It was complicated, so I broke it up

  ms = 2;
  ns = (1:ms-1);
  as = (-1).^(ns+1);
  bs = (as./(ns+1));
  
  cs = ((ms-1)./ns);

  ds = ns.*log2(ms);
  es = ns+1;
  fs = exp( -(ds./es)*EbNo );
  
  thrncoh = ((ms/2)/(ms-1)) * sum(bs.*((ms-1)./ns).*exp( -(ds./es)*EbNo ));

endfunction


function rxphi = fine_ex(timing_offset = 1,fsk_info)

  rand('state',1); 
  randn('state',1);

  bitcnt = 12051;
  test_bits = [zeros(1,100) rand(1,bitcnt)>.5]; %Random bits. Pad with zeros to prime the filters
  t_vec = [0 0 1 1];
  %test_bits = repmat(t_vec,1,ceil(24000/length(t_vec)));

  states = yafsk_init(fsk_info);
  Fs = states.Fs;
  Rb = states.Rs; 
  
  tx = yafsk_mod(states,test_bits);

  rx    = tx;
  rx    = rx(timing_offset:length(rx));

  [rx_bits states rxphi] = yafsk_demod_2a(states,rx);
  ber = 1;
  
  %thing to account for offset from input data to output data
  %No preamble detection yet
  ox = 1;
  for offset = (1:100)
    nerr = sum(xor(rx_bits(offset:length(rx_bits)),test_bits(1:length(rx_bits)+1-offset)));
    bern = nerr/(bitcnt-offset);
    if(bern < ber)
      ox = offset;
      best_nerr = nerr;
    end
    ber = min([ber bern]);
  end
  offset = ox;
  printf("\ncoarse timing: %d nerr: %d\n", offset, best_nerr);

endfunction

function phi=fine_2(aEsNodB,fsk_info,bits,offset,freq)
  [ber coh ncoh phi] = nfbert_2(aEsNodB,fsk_info, bits, offset, freq)
endfunction

function yafsk_rx_phi
  global fsk_setup
  pkg load parallel
  offrange = [1:2:101];

  setups = repmat(fsk_setup,1,length(offrange));
  phi = pararrayfun(nproc(),@fine_2,100,setups,1,offrange,0*ones(1,length(offrange)));
  
  close all;
  figure(1);
  clf;
  plot(offrange,real(phi),offrange,imag(phi));
  figure(2);
  plotyy(offrange,angle(phi),offrange,abs(phi))
endfunction

