% ofdm_lib.m
% David Rowe Mar 2017

#{
  Library of functions that implement a PSK OFDM modem.  Rate Fs
  verison of ofdm_rs.m with OFDM based up and down conversion, and all
  those nasty real-world details like fine freq, timing.  
#}

1;
qam16;

%-------------------------------------------------------------
% ofdm_init
%-------------------------------------------------------------

#{

  Modem frame has a pilot every Ns symbols. There are Ns-1 data
  symbols between every pilot.

   e.g. for Ns=4, Nc=6:
  
    |-Nc-|          Time
    DDDDDD           |
   PPPPPPPP  ---     |
    DDDDDD    |      |
    DDDDDD    Ns     |
    DDDDDD    |      |
   PPPPPPPP  ---    \|/ 
    DDDDDD    |      |

   Freq------------------>

   In this figure, time flows down, freq across.
#}

function states = ofdm_init(config)
  if isfield(config,"bps") bps = config.bps; else bps=2; end
  Rs = config.Rs; Tcp = config.Tcp; Ns = config.Ns; Nc = config.Nc;
  if isfield(config,"Np") Np = config.Np; else Np = 1; end
  if isfield(config,"Ntxtbits") Ntxtbits = config.Ntxtbits ; else Ntxtbits = 4; end
  if isfield(config,"Nuwbits") Nuwbits = config.Nuwbits ; else Nuwbits = 5*bps; end
  if isfield(config,"Nuwbits") Nuwbits = config.Nuwbits ; else Nuwbits = 5*bps; end
  if isfield(config,"ftwindow_width") ftwindow_width = config.ftwindow_width; else ftwindow_width = 11; end
  if isfield(config,"timing_mx_thresh") timing_mx_thresh = config.timing_mx_thresh; else timing_mx_thresh = 0.35; end
  if isfield(config,"tx_uw") tx_uw = config.tx_uw; else tx_uw = zeros(1,Nuwbits); end
  if isfield(config,"bad_uw_errors") bad_uw_errors = config.bad_uw_errors; else bad_uw_errors = 3; end
  
  states.Fs = 8000;
  states.bps = bps;
  states.Rs = Rs;
  states.Tcp = Tcp;
  states.Ns = Ns;                                 % one pilot every Ns symbols, e.g. Ns=3, ...PDDDPDDDP...
  states.Nc = Nc;                                 % Number of carriers
  states.M  = states.Fs/Rs;                       % oversampling rate
  states.Ncp = Tcp*states.Fs;
  states.Nbitsperframe = (Ns-1)*Nc*bps;           % total bits in all data symbols in modem frame
  states.Nsampersymbol = states.M+states.Ncp;     % number of samples in a single symbol
  states.Nsamperframe  = Ns*states.Nsampersymbol; % number of samples in a modem frame
  states.qam16 = [
    1 + j,  1 + j*3,  3 + j,  3 + j*3;
    1 - j,  1 - j*3,  3 - j,  3 - j*3;
   -1 + j, -1 + j*3, -3 + j, -3 + j*3;
   -1 - j, -1 - j*3, -3 - j, -3 - j*3]/3;
  states.Np = Np;                                 % number of modem frames per packet. In some modes we want
                                                  % the total packet of data to span multiple modem frames, e.g. HF data
                                                  % and/or when the FEC codeword is larger than the one
                                                  % modem frame.  In other modes (e.g. 700D/2020) Np=1, ie the modem frame
                                                  % is the same length as the packet/FEC frame.
  states.Nbitsperpacket = Np*states.Nbitsperframe;
  states.Tpacket = Np*Ns*(Tcp+1/Rs);              % time for one packet in ms

  states.Ntxtbits = Ntxtbits;                     % reserved bits/frame for auxillary text information.  Uncoded/unprotected so may
                                                  % be of limited use going forward, consider setting to 0
  states.Nuwbits  = Nuwbits;                      
  
  % some basic sanity checks
  assert(floor(states.M) == states.M);
  
  % UW symbol placement.  Use ofdm_dev.m, debug_false_sync() to test.
  % Note we need to fill each UW symbols with bits.  The LDPC decoder
  % works on symbols so we can't break up any symbols into UW/FEC
  % encoded bits.
  
  states.uw_ind = states.uw_ind_sym = [];
  for i=1:states.Nuwbits/bps
    ind_sym = floor(i*(Nc+1)/bps+1);
    states.uw_ind_sym = [states.uw_ind_sym ind_sym];   % symbol index
    for b=bps-1:-1:0
      states.uw_ind = [states.uw_ind bps*ind_sym-b];   % bit index
    end  
  end
  % how many of the first few frames have UW symbols in them
  Nsymsperframe = states.Nbitsperframe/states.bps;
  states.Nuwframes = ceil(states.uw_ind_sym(end)/Nsymsperframe);
  
  states.tx_uw = tx_uw;
  assert(length(states.tx_uw) == states.Nuwbits);
  tx_uw_syms = [];
  for b=1:bps:states.Nuwbits
    tx_uw_syms = [tx_uw_syms qpsk_mod(states.tx_uw(b:b+1))];
  end
  states.tx_uw_syms = tx_uw_syms;
  % if the UW has this many errors it is "bad", the binomal cdf can be used to set this:
  %   Nuw=12; plot(0:Nuw, binocdf(0:Nuw,Nuw,0.05)); hold on; plot(binocdf(0:Nuw,Nuw,0.5)); hold off;
  states.bad_uw_errors = bad_uw_errors;
  
  % use this to scale tx output to 16 bit short.  Adjusted by experiment
  % to have same RMS power as other FreeDV waveforms  
  states.amp_scale = 2E5*1.1491/1.06;

  % this is used to scale inputs to LDPC decoder to make it amplitude indep
  states.mean_amp = 0;

  % generate same pilots each time
  rand('seed',1);
  states.pilots = 1 - 2*(rand(1,Nc+2) > 0.5);
  %printf("number of pilots total: %d\n", length(states.pilots));
  
  % carrier tables for up and down conversion
  fcentre = 1500;
  alower = fcentre - Rs * (Nc/2);  % approx frequency of lowest carrier
  Nlower = round(alower / Rs) - 1; % round this to nearest integer multiple from 0Hz to keep DFT happy
  %printf("  fcentre: %f alower: %f alower/Rs: %f Nlower: %d\n", fcentre, alower, alower/Rs, Nlower);
  w = (Nlower:Nlower+Nc+1)*2*pi/(states.Fs/Rs);
  W = zeros(Nc+2,states.M);
  for c=1:Nc+2
    W(c,:) = exp(j*w(c)*(0:states.M-1));
  end
  states.w = w;
  states.W = W;

  % fine timing search +/- window_width/2 from current timing instant,
  % adjust this to be roughly the maximum delay spread
  states.ftwindow_width = ftwindow_width;

  % magic number we adjust by experiment (see ofdm_dev.m acquisition tests, blog post on 700D sync)
  states.timing_mx_thresh = timing_mx_thresh;

  % Receive buffer: D P DDD P DDD P DDD P D
  %                         ^
  % also see ofdm_demod() ...

  %                       D                 P DDD P DDD P DDD             P                    D
  states.Nrxbuf = states.Nsampersymbol + 3*states.Nsamperframe + states.Nsampersymbol + states.Nsampersymbol;
  states.rxbuf = zeros(1, states.Nrxbuf);
 
  % default settings on a bunch of options and states

  states.verbose = 0;
  states.timing_en = 1;
  states.foff_est_en = 1;
  states.phase_est_en = 1;
  states.phase_est_bandwidth = "high";
  states.dpsk = 0;
  
  states.foff_est_gain = 0.1;
  states.foff_est_hz = 0;
  states.sample_point = states.timing_est = 1;
  states.nin = states.Nsamperframe;
  states.timing_valid = 0;
  states.timing_mx = 0;
  states.coarse_foff_est_hz = 0;

  states.foff_metric = 0;
  
  % generate OFDM pilot symbol, used for timing and freq offset est

  rate_fs_pilot_samples = states.pilots * W/states.M;

  % During tuning it was found that not including the cyc prefix in
  % rate_fs_pilot_samples produced better fest results
  
  %states.rate_fs_pilot_samples = [rate_fs_pilot_samples(states.M-states.Ncp+1:states.M) rate_fs_pilot_samples];
  states.rate_fs_pilot_samples = [zeros(1,states.Ncp) rate_fs_pilot_samples];

  % pre-compute a constant used to detect valid modem frames

  Npsam = length(states.rate_fs_pilot_samples);
  states.timing_norm = Npsam*(states.rate_fs_pilot_samples * states.rate_fs_pilot_samples');
  % printf("timing_norm: %f\n", states.timing_norm)

  % sync state machine

  states.sync_state = states.last_sync_state = 'search'; 
  states.uw_errors = 0;
  states.sync_counter = 0;
  states.frame_count = 0;                                 % number of frames we have been in sync
  states.sync_start = 0;
  states.sync_end = 0;
  states.modem_frame = 0;                                 % keep track of how many frames received in packet
  
  % LDPC code is optionally enabled

  states.rate = 1.0;
  states.ldpc_en = 0;

  % init some output states for logging
  
  states.rx_sym = zeros(1+Ns+1+1, Nc+2);

  % Es/No (SNR) est states
  
  states.noise_var = 0;
  states.sig_var = 0;
  states.clock_offset_est = 0;

  % automated tests
  test_qam16_mod_demod(states.qam16);
  test_assemble_disassemble(states);
endfunction


%------------------------------------------------------------------------------
% ofdm_init_mode - Helper function to set up modems for various FreeDV modes,
%                  and parse mode string.
%------------------------------------------------------------------------------

function config = ofdm_init_mode(mode="700D")
  Tcp = 0.002; Ns=8;

  % some "canned" modes
  if strcmp(mode,"700D")
    Ts = 0.018; Nc = 17;
  elseif strcmp(mode,"2020")
    Ts = 0.0205; Nc = 31;
  elseif strcmp(mode,"2200")
    Tframe = 0.175; Ts = Tframe/Ns; Nc = 37;
  elseif strcmp(mode,"qam16")
    Ns=5; config.Np=5; Tcp = 0.004; Ts = 0.016; Nc = 33;
    config.bps=4; config.Ntxtbits = 0; config.Nuwbits = 15*4; config.bad_uw_errors = 5;
    config.ftwindow_width = 32;
  elseif strcmp(mode,"datac1")
    Ns=5; config.Np=18; Tcp = 0.006; Ts = 0.016; Nc = 18;
    config.Ntxtbits = 0; config.Nuwbits = 12; config.bad_uw_errors = 2;
    config.ftwindow_width = 32;
  elseif strcmp(mode,"datac2")
    Ns=5; config.Np=36; Tcp = 0.006; Ts = 0.016; Nc = 9;
    config.Ntxtbits = 0; config.Nuwbits = 12; config.bad_uw_errors = 1;
    config.ftwindow_width = 32;
  elseif strcmp(mode,"datac3")
    Ns=5; config.Np=11; Tcp = 0.006; Ts = 0.016; Nc = 9;
    config.Ntxtbits = 0; config.Nuwbits = 24; config.bad_uw_errors = 5;
    config.ftwindow_width = 32; config.timing_mx_thresh = 0.30;
    config.tx_uw = [1 1 0 0  1 0 1 0  1 1 1 1  0 0 0 0  1 1 1 1  0 0 0 0];
  elseif strcmp(mode,"1")
    Ns=5; config.Np=10; Tcp=0; Tframe = 0.1; Ts = Tframe/Ns; Nc = 1;
  else
    % try to parse mode string for user defined mode
    vec = sscanf(mode, "Ts=%f Nc=%d Ncp=%f");
    Ts=vec(1); Nc=vec(2); Ncp=vec(3);
  end
  Rs=1/Ts;
  config.Rs = Rs; config.Tcp = Tcp; config.Ns = Ns; config.Nc = Nc;
end


%------------------------------------------------------------------------------
% print_config - utility function to use ascsii-art to describe the modem frame
%------------------------------------------------------------------------------

function print_config(states)
  ofdm_load_const;
  printf("Nc=%d Ts=%4.3f Tcp=%4.3f Ns: %d Np: %d\n", Nc, 1/Rs, Tcp, Ns, Np);
  printf("Nsymperframe: %d Nbitsperpacket: %d Nsamperframe: %d Ntxtbits: %d Nuwbits: %d Nuwframes: %d\n",
          Ns*Nc, Nbitsperpacket, Nsamperframe, Ntxtbits, Nuwbits, Nuwframes);
  printf("uncoded bits/s: %4.1f\n",  Nbitsperpacket*Fs/(Np*Nsamperframe));
  s=1; u=1; Nuwsyms=length(uw_ind_sym);
  for f=1:Np
    for r=1:Ns
      for c=1:Nc+2
        if r == 1
          sym="P";
        elseif c>1 && c <=(Nc+1)
          sym=".";
          if (u <= Nuwsyms) && (s == uw_ind_sym(u)) sym="U"; u++; end
          s++;
        else
          sym=" ";
        end
        printf("%s",sym);
      end
      printf("\n");
    end
  end  
end

% Gray coded QPSK modulation function
function symbol = qpsk_mod(two_bits)
    two_bits_decimal = sum(two_bits .* [2 1]); 
    switch(two_bits_decimal)
        case (0) symbol =  1;
        case (1) symbol =  j;
        case (2) symbol = -j;
        case (3) symbol = -1;
    endswitch
endfunction

% Gray coded QPSK demodulation function
function two_bits = qpsk_demod(symbol)
    bit0 = real(symbol*exp(j*pi/4)) < 0;
    bit1 = imag(symbol*exp(j*pi/4)) < 0;
    two_bits = [bit1 bit0];
endfunction


function out = freq_shift(in, foff, Fs)
  foff_rect = exp(j*2*pi*foff/Fs);
  foff_phase_rect = exp(j*0);
  
  for r=1:length(in)
    foff_phase_rect *= foff_rect;
    out(r) = in(r)*foff_phase_rect;
  end
endfunction


% -----------------------------------------------------------------
% ofdm_mod - modulates a complete packet (one or more modem frames)
% ----------------------------------------------------------------

function tx = ofdm_mod(states, tx_bits)
  ofdm_load_const;
  assert(length(tx_bits) == Nbitsperpacket);
 
  % map to symbols in linear array

  if bps == 1
    tx_sym_lin = 2*tx_bits - 1;
  end
  if bps == 2
    for s=1:Nbitsperpacket/bps
      tx_sym_lin(s) = qpsk_mod(tx_bits(2*(s-1)+1:2*s));
    end
  end  
  if bps == 4
    for s=1:Nbitsperpacket/bps
      tx_sym_lin(s) = qam16_mod(states.qam16,tx_bits(4*(s-1)+1:4*s))*exp(-j*pi/4);
    end
  end
  
  tx = ofdm_txframe(states, tx_sym_lin);
endfunction


% ----------------------------------------------
% ofdm_txframe - modulates one packet of symbols
% ----------------------------------------------

function tx = ofdm_txframe(states, tx_sym_lin)
  ofdm_load_const;
  assert(length(tx_sym_lin) == Nbitsperpacket/bps);

  % place data symbols in multi-carrier frame with pilots and boundary carriers

  s = 1; tx_frame = zeros(Np*Ns,Nc+2);
  for r=1:Np*Ns
    if mod(r-1,Ns) == 0
      % row of pilots
      tx_frame(r,:) = pilots;
    else
      % row of data symbols
      arowofsymbols = tx_sym_lin(s:s+Nc-1);
      tx_frame(r,2:Nc+1) = arowofsymbols;
      s += Nc;
      if states.dpsk
        tx_frame(r,2:Nc+1) = tx_frame(r,2:Nc+1) .* tx_frame(r-1,2:Nc+1);
      end
    end
  end
  % make sure we use all the symbols
  assert((s-1) == length(tx_sym_lin));

  % OFDM upconvert symbol by symbol so we can add CP

  tx = [];
  for r=1:Ns*Np
    asymbol = tx_frame(r,:) * W/M;
    asymbol_cp = [asymbol(M-Ncp+1:M) asymbol];
    tx = [tx asymbol_cp];
  end
endfunction


% -----------------------------------------------------------
% est_timing
% -----------------------------------------------------------

#{
  Correlates the OFDM pilot symbol samples with a window of received
  samples to determine the most likely timing offset.  Combines two
  frames pilots so we need at least Nsamperframe+M+Ncp samples in rx.

  Can be used for acquisition (coarse timing), and fine timing.  Tends
  to break down when freq offset approaches +/- symbol rate (e.g +/-
  25 Hz for 700D).
#}

function [t_est timing_valid timing_mx av_level] = est_timing(states, rx, rate_fs_pilot_samples, step)
    ofdm_load_const;
    Npsam = length(rate_fs_pilot_samples);
    
    Ncorr = length(rx) - (Nsamperframe+Npsam);
    corr = zeros(1,Ncorr);
    %printf("Npsam: %d M+Ncp: %d Ncorr: %d Nsamperframe: %d step: %d\n", Npsam,  M+Ncp, Ncorr, Nsamperframe, step);
    
    % normalise correlation so we can compare to a threshold across varying input levels

    av_level = 2*sqrt(states.timing_norm*(rx*rx')/length(rx)) + 1E-12;

    % correlate with pilots at start and end of frame to determine timing offset
    
    for i=1:step:Ncorr
      rx1     = rx(i:i+Npsam-1); rx2 = rx(i+Nsamperframe:i+Nsamperframe+Npsam-1);
      corr_st = rx1 * rate_fs_pilot_samples'; corr_en = rx2 * rate_fs_pilot_samples';
      corr(i) = (abs(corr_st) + abs(corr_en))/av_level;
    end

    [timing_mx t_est] = max(corr);
    % only declare timing valid if there are enough samples in rxbuf to demodulate a frame
    timing_valid = (abs(rx(t_est)) > 0) && (timing_mx > timing_mx_thresh);
    
    if verbose > 1
      printf("  av_level: %5.4f mx: %4.3f timing_est: %4d timing_valid: %d\n", av_level, timing_mx, t_est, timing_valid);
    end
    if verbose > 2
      figure(3); clf;
      subplot(211); plot(rx)
      subplot(212); plot(corr)
      figure(4); clf; plot(real(rate_fs_pilot_samples));
    end

endfunction


% -----------------------------------------------------------
% est_freq_offset
% -----------------------------------------------------------

#{
  Determines frequency offset at current timing estimate, used for
  coarse freq offset estimation during acquisition.

  This estimator works well for AWGN channels but has problems with
  fading channels.  With stationary/slow fading channels (say a notch
  in the spectrum), ot exhibits bias which can delay sync for 10's of
  seconds.
#}

function [foff_est states] = est_freq_offset(states, rx, rate_fs_pilot_samples, t_est)
    ofdm_load_const;
    Npsam = length(rate_fs_pilot_samples);

    % Freq offset can be considered as change in phase over two halves
    % of pilot symbols.  We average this statistic over this and next
    % frames pilots.

    Npsam2 = floor(Npsam/2);
    p1 = rx(t_est:t_est+Npsam2-1) * rate_fs_pilot_samples(1:Npsam2)';
    p2 = rx(t_est+Npsam2:t_est+Npsam-1) * rate_fs_pilot_samples(Npsam2+1:Npsam)';
    p3 = rx(t_est+Nsamperframe:t_est+Nsamperframe+Npsam2-1) * rate_fs_pilot_samples(1:Npsam2)';
    p4 = rx(t_est+Nsamperframe+Npsam2:t_est+Nsamperframe+Npsam-1) * rate_fs_pilot_samples(Npsam2+1:Npsam)';
   
    Fs1 = Fs/(Npsam/2);

    states.foff_metric = (conj(p1)*p2 + conj(p3)*p4);
    foff_est = Fs1*angle(states.foff_metric)/(2*pi);

    if states.verbose > 1
        printf("  foff_metric: %f %f foff_est: %f\n", real(states.foff_metric), imag(states.foff_metric), foff_est);
    end
 
endfunction


% -----------------------------------------------------------
% est_freq_offset_pilot_corr
% -----------------------------------------------------------

#{
  Determines frequency offset at current timing estimate, used for
  coarse freq offset estimation during acquisition.

  This is an alternative algorithm to est_freq_offset() above that is less noisey
  and performs better on HF channels using the acquistion tests in ofdm_dev.m
#}

function foff_est = est_freq_offset_pilot_corr(states, rx, rate_fs_pilot_samples, t_est)
    ofdm_load_const;
    Npsam = length(rate_fs_pilot_samples);

    % extract pilot samples from either end of frame
    rx1  = rx(t_est:t_est+Npsam-1); rx2 = rx(t_est+Nsamperframe:t_est+Nsamperframe+Npsam-1);

    % "mix" these down (correlate) with 0 Hz offset pilot samples
    corr_st = rx1 .* conj(rate_fs_pilot_samples);
    corr_en = rx2 .* conj(rate_fs_pilot_samples);

    % sample sum of DFT magnitude of correlated signals at each freq offset and look for peak
    st = -20; en = 20; foff_est = 0; Cabs_max = 0;

    for f=st:en
       w = 2*pi*f/Fs;
       C_st = corr_st * exp(j*w*(0:Npsam-1))';
       C_en = corr_en * exp(j*w*(0:Npsam-1))';
       Cabs = abs(C_st) + abs(C_en);
       if Cabs > Cabs_max
         Cabs_max = Cabs;
         foff_est = f;
       end
    end
    
    if states.verbose > 1
      printf("  foff_est: %f\n", foff_est);
    end
    if verbose > 2
      figure(10); clf;
      plot(st:en,C(Fs/2+st:Fs/2+en)); grid;
    end 
endfunction


% ----------------------------------------------------------------------------------
% ofdm_sync_search - attempts to find coarse sync parameters for modem initial sync
% ----------------------------------------------------------------------------------

function [timing_valid states] = ofdm_sync_search(states, rxbuf_in)
  ofdm_load_const;

  % insert latest input samples into rxbuf so it is primed for when we have to call ofdm_demod()

  states.rxbuf(1:Nrxbuf-states.nin) = states.rxbuf(states.nin+1:Nrxbuf);
  states.rxbuf(Nrxbuf-states.nin+1:Nrxbuf) = rxbuf_in;

  % Attempt coarse timing estimate (i.e. detect start of frame) at a range of frequency offsets

  st = M+Ncp + Nsamperframe + 1; en = st + 2*Nsamperframe +  M+Ncp - 1;
  timing_mx = 0; fcoarse = 0; timing_valid = 0; ct_est = 1;
  for afcoarse=-40:40:40
    % vector of local oscillator samples to shift input vector
    % these could be computed on the fly to save memory, or pre-computed in flash at tables as they are static

    if afcoarse != 0
      w = 2*pi*afcoarse/Fs;
      wvec = exp(-j*w*(0:2*Nsamperframe+M+Ncp-1));

      % choose best timing offset metric at this freq offset
      [act_est atiming_valid atiming_mx] = est_timing(states, wvec .* states.rxbuf(st:en), states.rate_fs_pilot_samples, 2);
    else
      % exp(-j*0) is just 1 when afcoarse is 0
      [act_est atiming_valid atiming_mx] = est_timing(states, states.rxbuf(st:en), states.rate_fs_pilot_samples, 2);
    end
    
    %printf("afcoarse: %f atiming_mx: %f\n", afcoarse, atiming_mx);
    
    if atiming_mx > timing_mx
      ct_est = act_est;
      timing_valid = atiming_valid;
      timing_mx = atiming_mx;
      fcoarse = afcoarse;
    end
  end
  
  % refine freq est within -/+ 20 Hz window  

  if fcoarse != 0
    w = 2*pi*fcoarse/Fs;
    wvec = exp(-j*w*(0:2*Nsamperframe+M+Ncp-1));
    foff_est = est_freq_offset_pilot_corr(states, wvec .* states.rxbuf(st:en), states.rate_fs_pilot_samples, ct_est);
    foff_est += fcoarse;
  else
    % exp(-j*0) is just 1 when fcoarse is 0
    foff_est = est_freq_offset_pilot_corr(states, states.rxbuf(st:en), states.rate_fs_pilot_samples, ct_est);
  end
 
  if verbose
    printf(" ct_est: %4d mx: %3.2f coarse_foff: %5.1f timing_valid: %d", ct_est, timing_mx, foff_est, timing_valid);
  end

  if timing_valid
    % potential candidate found ....

    % calculate number of samples we need on next buffer to get into sync

    states.nin = ct_est - 1;

    % reset modem states

    states.sample_point = states.timing_est = 1;
    states.foff_est_hz = foff_est;
  else
    states.nin = Nsamperframe;
  end
  
  states.timing_valid = timing_valid;
  states.timing_mx = timing_mx;
  states.coarse_foff_est_hz = foff_est;
endfunction


% ------------------------------------------
% ofdm_demod - Demodulates one frame of bits
% ------------------------------------------

#{ 

  For phase estimation we need to maintain buffer of 3 frames plus
  one pilot, so we have 4 pilots total. '^' is the start of current
  frame that we are demodulating.
           
  P DDD P DDD P DDD P
        ^
    
  Then add one symbol either side to account for movement in
  sampling instant due to sample clock differences:

  D P DDD P DDD P DDD P D
          ^

  Returns:
    rx_bits    - (hard decoded/raw/uncoded) demodulated data bits from packet
    aphase_est - phase est for each data symbol
    rx_np      - output data symbols after phase correction
    rx_amp     - amplitude estimates for each symbol
#}

function [states rx_bits aphase_est_pilot_log rx_np rx_amp] = ofdm_demod(states, rxbuf_in)
  ofdm_load_const;

  % insert latest input samples into rxbuf

  rxbuf(1:Nrxbuf-states.nin) = rxbuf(states.nin+1:Nrxbuf);
  rxbuf(Nrxbuf-states.nin+1:Nrxbuf) = rxbuf_in;

  % get latest freq offset estimate

  woff_est = 2*pi*foff_est_hz/Fs;

  % update timing estimate --------------------------------------------------

  delta_t = coarse_foff_est_hz = timing_valid = timing_mx = 0;
  if timing_en
    % update timing at start of every frame

    st = M+Ncp + Nsamperframe + 1 - floor(ftwindow_width/2) + (timing_est-1);
    en = st + Nsamperframe-1 + M+Ncp + ftwindow_width-1;
          
    [ft_est timing_valid timing_mx] = est_timing(states, rxbuf(st:en) .* exp(-j*woff_est*(st:en)), rate_fs_pilot_samples, 1);
    % printf("  timing_est: %d ft_est: %d timing_valid: %d timing_mx: %d\n", timing_est, ft_est, timing_valid, timing_mx);
    
    if timing_valid
      timing_est = timing_est + ft_est - ceil(ftwindow_width/2);

      % Black magic to keep sample_point inside cyclic prefix.  Or something like that.

      delta_t = ft_est - ceil(ftwindow_width/2);
      sample_point = max(timing_est+Ncp/4, sample_point);
      sample_point = min(timing_est+Ncp, sample_point);
    end
    
    if verbose > 1
      printf("  ft_est: %2d mx: %3.2f coarse_foff: %4.1f foff: %4.1f\n", ft_est, timing_mx, coarse_foff_est_hz, foff_est_hz);
    end

  end

  % down convert at current timing instant----------------------------------

  % todo: this cld be more efficent, as pilot r becomes r-Ns on next frame

  rx_sym = zeros(1+Ns+1+1, Nc+2);

  % previous pilot
  
  st = M+Ncp + Nsamperframe + (-Ns)*(M+Ncp) + 1 + sample_point; en = st + M - 1;

  for c=1:Nc+2
    acarrier = rxbuf(st:en) .* exp(-j*woff_est*(st:en)) .* conj(W(c,:));
    rx_sym(1,c) = sum(acarrier);
  end

  % pilot - this frame - pilot

  for rr=1:Ns+1 
    st = M+Ncp + Nsamperframe + (rr-1)*(M+Ncp) + 1 + sample_point; en = st + M - 1;
    for c=1:Nc+2
      acarrier = rxbuf(st:en) .* exp(-j*woff_est*(st:en)) .* conj(W(c,:));
      rx_sym(rr+1,c) = sum(acarrier);
    end
  end

  % next pilot

  st = M+Ncp + Nsamperframe + (2*Ns)*(M+Ncp) + 1 + sample_point; en = st + M - 1;
  for c=1:Nc+2
    acarrier = rxbuf(st:en) .* exp(-j*woff_est*(st:en)) .* conj(W(c,:));
    rx_sym(Ns+3,c) = sum(acarrier);
  end
      
  % est freq err based on all carriers ------------------------------------
      
  if foff_est_en
    freq_err_rect = sum(rx_sym(2,:))' * sum(rx_sym(2+Ns,:));

    % prevent instability in atan(im/re) when real part near 0 

    freq_err_rect += 1E-6;

    %printf("freq_err_rect: %f %f angle: %f\n", real(freq_err_rect), imag(freq_err_rect), angle(freq_err_rect));
    freq_err_hz = angle(freq_err_rect)*Rs/(2*pi*Ns);
    foff_est_hz = foff_est_hz + foff_est_gain*freq_err_hz;
  end

  % OK - now channel for each carrier and correct phase  ----------------------------------

  achannel_est_rect = zeros(1,Nc+2);
  for c=2:Nc+1

    % estimate channel for this carrier using an average of 12 pilots
    % in a rect 2D window centred on this carrier
    
    % PPP  <-- frame-1
    % ---
    % PPP  <-- you are here
    % DDD
    % DDD
    % PPP  <-- frame+1
    % ---
    % PPP  <-- frame+2
    
    if isfield(states, "phase_est_bandwidth")
      phase_est_bandwidth = states.phase_est_bandwidth;
    else
      phase_est_bandwidth = "low";
    end
    
    if strcmp(phase_est_bandwidth, "high")
      % Only use pilots at start and end of this frame to track quickly changes in phase
      % present.  Useful for initial sync where freq offset est may be a bit off, and
      % for high Doppler channels.  As less pilots are averaged, low SNR performance
      % will be poorer.
      achannel_est_rect(c) =  sum(rx_sym(2,c)*pilots(c)');      % frame    
      achannel_est_rect(c) += sum(rx_sym(2+Ns,c)*pilots(c)');   % frame+1
    else
      % Average over a bunch of pilots in adjacent carriers, and past and future frames, good
      % low SNR performance, but will fall over with high Doppler of freq offset.
      cr = c-1:c+1;
      achannel_est_rect(c) =  sum(rx_sym(2,cr)*pilots(cr)');      % frame    
      achannel_est_rect(c) += sum(rx_sym(2+Ns,cr)*pilots(cr)');   % frame+1

      % use next step of pilots in past and future

      achannel_est_rect(c) += sum(rx_sym(1,cr)*pilots(cr)');      % frame-1
      achannel_est_rect(c) += sum(rx_sym(2+Ns+1,cr)*pilots(cr)'); % frame+2
    end
  end
  
  if strcmp(phase_est_bandwidth, "high")
    achannel_est_rect /= 2;
  else
    achannel_est_rect /= 12;
  end
  
  % pilots are estimated over 12 pilot symbols, so find average

  aphase_est_pilot = angle(achannel_est_rect);
  aamp_est_pilot = abs(achannel_est_rect);

  % correct phase offset using phase estimate, and demodulate
  % bits, separate loop as it runs across cols (carriers) to get
  % frame bit ordering correct

  aphase_est_pilot_log = [];
  rx_bits = []; rx_np = []; rx_amp = [];
  for rr=1:Ns-1
    for c=2:Nc+1
      if phase_est_en
        if states.dpsk
          rx_corr = rx_sym(rr+2,c) *  rx_sym(rr+1,c)';
        else
          rx_corr = rx_sym(rr+2,c) * exp(-j*aphase_est_pilot(c));
        end
      else
        rx_corr = rx_sym(rr+2,c);
      end
      rx_np = [rx_np rx_corr];
      rx_amp = [rx_amp aamp_est_pilot(c)];
      if bps == 1
        abit = real(rx_corr) > 0;
      end
      if bps == 2
        abit = qpsk_demod(rx_corr);
      end
      if bps == 4
        abit = qam16_demod(states.qam16, rx_corr*exp(j*pi/4));
      end
      rx_bits = [rx_bits abit];
    end % c=2:Nc+1
    aphase_est_pilot_log = [aphase_est_pilot_log; aphase_est_pilot(2:Nc+1)];
  end 

  % Adjust nin to take care of sample clock offset

  nin = Nsamperframe;
  if timing_en && timing_valid
    states.clock_offset_est = 0.9*states.clock_offset_est + 0.1*abs(states.timing_est - timing_est)/Nsamperframe;
    thresh = (M+Ncp)/8;
    tshift = (M+Ncp)/4;
    if timing_est > thresh
      nin = Nsamperframe+tshift;
      timing_est -= tshift;
      sample_point -= tshift;
    end
    if timing_est < -thresh
      nin = Nsamperframe-tshift;
      timing_est += tshift;
      sample_point += tshift;
    end
  end

  % estimates of signal and noise power (see cohpsk.m for further explanation)
  % signal power is distance from axis on complex plane
  % we just measure noise power on imag axis, as it isn't affected by fading
  % using all symbols in frame worked better than just pilots
  
  sig_var = sum(abs(rx_np) .^ 2)/length(rx_np);
  sig_rms = sqrt(sig_var);
  
  sum_x = 0;
  sum_xx = 0;
  n = 0;
  for i=1:length(rx_np)
    s = rx_np(i);
    if abs(real(s)) > sig_rms 
      % select two constellation points on real axis
      sum_x  += imag(s);
      sum_xx += imag(s)*imag(s);
      n++;
    end
  end
   
  noise_var = 0;
  if n > 1
    noise_var = (n*sum_xx - sum_x*sum_x)/(n*(n-1));
  end

  % Total noise power is twice estimate of imaginary-axis noise.  This
  % effectively gives us the an estimate of Es/No
  
  states.noise_var = 2*noise_var; 
  states.sig_var = sig_var;

  % maintain mean amp estimate for LDPC decoder

  states.mean_amp = 0.9*states.mean_amp + 0.1*mean(rx_amp);

  states.achannel_est_rect = achannel_est_rect;
  states.rx_sym = rx_sym;
  states.rxbuf = rxbuf;
  states.nin = nin;
  states.timing_valid = timing_valid;
  states.timing_mx = timing_mx;
  states.timing_est = timing_est;
  states.sample_point = sample_point;
  states.delta_t = delta_t;
  states.foff_est_hz = foff_est_hz;
  states.coarse_foff_est_hz = coarse_foff_est_hz; % just used for tofdm
endfunction


% ----------------------------------------------------------------------------------
% assemble_modem_packet - assemble modem packet from UW, payload, and txt bits
% ----------------------------------------------------------------------------------

function modem_frame = assemble_modem_packet(states, payload_bits, txt_bits)
  ofdm_load_const;

  # Due to the operation of the FEC encoder or interleaver, Tx data
  # usually comes in "packet size" chunks, so assembly operates on an
  # entire packet (multiple modem frames if Np>1)

  p = 1; u = 1;
  modem_frame = zeros(1,Nbitsperpacket);
  
  for b=1:Nbitsperpacket-Ntxtbits;
    if (u <= Nuwbits) && (b == uw_ind(u))
      modem_frame(b) = tx_uw(u++);
    else
      modem_frame(b) = payload_bits(p++);
    end  
  end
  t = 1;
  for b=Nbitsperpacket-Ntxtbits+1:Nbitsperpacket
    modem_frame(b) = txt_bits(t++);
  end
  assert(u == (Nuwbits+1));
  assert(p = (length(payload_bits)+1));
endfunction


% ----------------------------------------------------------------------------------
% assemble_modem_packet_symbols - assemble modem packet from UW, payload, and txt bits
% ----------------------------------------------------------------------------------

function modem_frame = assemble_modem_packet_symbols(states, payload_syms, txt_syms)
  ofdm_load_const;

  Nsymsperpacket = Nbitsperpacket/bps;
  Nuwsyms = Nuwbits/bps;
  Ntxtsyms = Ntxtbits/bps;
  modem_frame = zeros(1,Nsymsperpacket);
  p = 1; u = 1;

  for s=1:Nsymsperpacket-Ntxtsyms;
    if (u <= Nuwsyms) && (s == uw_ind_sym(u))
      modem_frame(s) = states.tx_uw_syms(u++);
    else
      modem_frame(s) = payload_syms(p++);
    end
  end
  t = 1;
  for s=Nsymsperpacket-Ntxtsyms+1:Nsymsperpacket
    modem_frame(s) = txt_syms(t++);
  end
  assert(u == (Nuwsyms+1));
  assert(p = (length(payload_syms)+1));
endfunction


% ------------------------------------------------------------------------------------------------
% extract_uw - extract just the UW from the first few frames of a packet, to check UW
%              during acquisition
% -------------------------------------------------------------------------------------------------

function rx_uw = extract_uw(states, rx_syms)
  ofdm_load_const;

  Nsymsperframe = Nbitsperframe/bps;
  assert(length(rx_syms) == Nuwframes*Nsymsperframe);
  Nuwsyms = Nuwbits/bps;
  rx_uw_syms = zeros(1,Nuwsyms);
  u = 1;
 
  for s=1:Nuwframes*Nsymsperframe
    if (u <= Nuwsyms) && (s == uw_ind_sym(u))
      rx_uw_syms(u++) = rx_syms(s);
    end
  end
  assert(u == (Nuwsyms+1));

  % now demodulate UW bits  
  rx_uw = zeros(1,Nuwbits);
  
  for s=1:Nuwsyms
    if bps == 2
      rx_uw(bps*(s-1)+1:bps*s) = qpsk_demod(rx_uw_syms(s));
    elseif bps == 4
      rx_uw(bps*(s-1)+1:bps*s) = qam16_demod(states.qam16,rx_uw_syms(s)*exp(j*pi/4));
    end
  end
endfunction


% ------------------------------------------------------------------------------------------------
% disassemble_modem_packet - extract UW, txt bits, and payload symbols from a packet of symbols
% -------------------------------------------------------------------------------------------------

function [rx_uw payload_syms payload_amps txt_bits] = disassemble_modem_packet(states, modem_frame_syms, modem_frame_amps)
  ofdm_load_const;

  Nsymsperpacket = Nbitsperpacket/bps;
  Nuwsyms = Nuwbits/bps;
  Ntxtsyms = Ntxtbits/bps;
  payload_syms = zeros(1,Nsymsperpacket-Nuwsyms-Ntxtsyms);
  payload_amps = zeros(1,Nsymsperpacket-Nuwsyms-Ntxtsyms);
  rx_uw_syms = zeros(1,Nuwsyms);
  txt_syms = zeros(1,Ntxtsyms);
  p = 1; u = 1;
 
  for s=1:Nsymsperpacket-Ntxtsyms;
    if (u <= Nuwsyms) && (s == uw_ind_sym(u))
      rx_uw_syms(u++) = modem_frame_syms(s);
    else
      payload_syms(p) = modem_frame_syms(s);
      payload_amps(p++) = modem_frame_amps(s);
    end
  end
  t = 1;
  for s=Nsymsperpacket-Ntxtsyms+1:Nsymsperpacket
    txt_syms(t++) = modem_frame_syms(s);
  end
  assert(u == (Nuwsyms+1));
  assert(p = (Nsymsperpacket+1));

  % now demodulate UW and txt bits
  
  rx_uw = zeros(1,Nuwbits);
  txt_bits = zeros(1,Ntxtbits);
  
  for s=1:Nuwsyms
    if bps == 2
      rx_uw(bps*(s-1)+1:bps*s) = qpsk_demod(rx_uw_syms(s));
    elseif bps == 4
      rx_uw(bps*(s-1)+1:bps*s) = qam16_demod(states.qam16,rx_uw_syms(s)*exp(j*pi/4));
    end
  end
  for s=1:Ntxtsyms
    txt_bits(2*s-1:2*s) = qpsk_demod(txt_syms(s));
  end

endfunction


%-----------------------------------------------------------------------
% ofdm_rand - a psuedo-random number generator that we can implement
%             in C with identical results to Octave.  Returns an unsigned
%             int between 0 and 32767
%-----------------------------------------------------------------------

function r = ofdm_rand(n)
  r = zeros(1,n); seed = 1;
  for i=1:n
    seed = mod(1103515245 * seed + 12345, 32768);
    r(i) = seed;
  end
endfunction


%-----------------------------------------------------------------------
% create_ldpc_test_frame - generate a test frame of bits
%-----------------------------------------------------------------------

function [tx_bits payload_data_bits codeword] = create_ldpc_test_frame(states, coded_frame=1)
  ofdm_load_const;
  ldpc;
  gp_interleaver;
  
  if coded_frame
    % Set up LDPC code

    mod_order = 4; bps = 2; modulation = 'QPSK'; mapping = 'gray';

    init_cml('~/cml/'); % TODO: make this path sensible and portable
    load HRA_112_112.txt
    [code_param framesize rate] = ldpc_init_user(HRA_112_112, modulation, mod_order, mapping);
    assert(Nbitsperframe == (code_param.coded_bits_per_frame + Nuwbits + Ntxtbits));

    payload_data_bits = round(ofdm_rand(code_param.data_bits_per_frame)/32767);
    codeword = LdpcEncode(payload_data_bits, code_param.H_rows, code_param.P_matrix);
    Nsymbolsperframe = length(codeword)/bps;
  
    % need all these steps to get actual raw codeword bits at demod ..
  
    tx_symbols = [];
    for s=1:Nsymbolsperframe
      tx_symbols = [tx_symbols qpsk_mod( codeword(2*(s-1)+1:2*s) )];
    end

    tx_symbols = gp_interleave(tx_symbols);

    codeword_raw = [];
    for s=1:Nsymbolsperframe
      codeword_raw = [codeword_raw qpsk_demod(tx_symbols(s))];
    end
  else
    codeword_raw = round(ofdm_rand(Nbitsperpacket-(Nuwbits+Ntxtbits))/32767);
  end
  
  % insert UW and txt bits
  
  tx_bits = assemble_modem_packet(states, codeword_raw, zeros(1,Ntxtbits));
  assert(Nbitsperpacket == length(tx_bits));

endfunction

% automated test

function test_assemble_disassemble(states)
  ofdm_load_const;

  Nsymsperpacket = Nbitsperpacket/bps;
  Ndatabitsperpacket = Nbitsperpacket-(Nuwbits+Ntxtbits);
  Ndatasymsperpacket = Ndatabitsperpacket/bps;
  codeword_bits = round(ofdm_rand(Ndatabitsperpacket)/32767);
  tx_bits = assemble_modem_packet(states, codeword_bits, zeros(1,Ntxtbits));

  tx_syms = zeros(1,Nsymsperpacket);
  for s=1:Nsymsperpacket
    if bps == 2
      tx_syms(s) = qpsk_mod(tx_bits(bps*(s-1)+1:bps*s));
    elseif bps == 4
      tx_syms(s) = qam16_mod(states.qam16,tx_bits(bps*(s-1)+1:bps*s));
    end
  end
  codeword_syms = zeros(1,Ndatasymsperpacket);
  for s=1:Ndatasymsperpacket
    if bps == 2
      codeword_syms(s) = qpsk_mod(codeword_bits(bps*(s-1)+1:bps*s));
    elseif bps == 4
      codeword_syms(s) = qam16_mod(states.qam16,codeword_bits(bps*(s-1)+1:bps*s));
    end
  end

  [rx_uw rx_codeword_syms payload_amps txt_bits] = disassemble_modem_packet(states, tx_syms, ones(1,Nsymsperpacket));
  assert(rx_uw == states.tx_uw);
  Ndatasymsperframe = (Nbitsperpacket-(Nuwbits+Ntxtbits))/bps;
  assert(codeword_syms == rx_codeword_syms);
endfunction

%-------------------------------------------------------
% sync_state_machine - determines sync state based on UW
%                      700D/2020 version
%-------------------------------------------------------

#{
  Due to the low pilot symbol insertion rate and acquisition issues
  the earlier OFDM modem waveforms (700D and 2020) need a complex
  state machine to help them avoid false sync.
#}

function states = sync_state_machine(states, rx_uw)
  ofdm_load_const;
  next_state = states.sync_state;
  states.sync_start = states.sync_end = 0;
  
  if strcmp(states.sync_state,'search') 

    if states.timing_valid
      states.frame_count = 0;
      states.sync_counter = 0;
      states.modem_frame = 0;
      states.sync_start = 1;
      next_state = 'trial';
    end
  end
        
  if strcmp(states.sync_state,'synced') || strcmp(states.sync_state,'trial')

    states.frame_count++;

    % UW occurs at the start of a packet
    if states.modem_frame == 0
        states.uw_errors = sum(xor(tx_uw,rx_uw));

        if strcmp(states.sync_state,'trial')
          if states.uw_errors >= states.bad_uw_errors
            states.sync_counter++;
            states.frame_count = 0;
          end
          if states.sync_counter == 2
            next_state = "search";
            states.phase_est_bandwidth = "high";
          end
          if states.frame_count == 4
            next_state = "synced";
            % change to low bandwidth, but more accurate phase estimation
            states.phase_est_bandwidth = "low";
          end
          if states.uw_errors < 2
            next_state = "synced";
            % change to low bandwidth, but more accurate phase estimation
            states.phase_est_bandwidth = "low";
          else
            next_state = "search"
          end  
        end

        if strcmp(states.sync_state,'synced')
          if states.uw_errors > 2
            states.sync_counter++;
          else
            states.sync_counter = 0;
          end

          if states.sync_counter == 6
            next_state = "search";
            states.phase_est_bandwidth = "high";
          end
        end
      end % if modem_frame == 0 ....

      % keep track of where we are up to in packet
      states.modem_frame++;
      if (states.modem_frame >= states.Np) states.modem_frame = 0; end
  end
  
  states.last_sync_state = states.sync_state;
  states.sync_state = next_state;
endfunction


%-------------------------------------------------------
% sync_state_machine_data - data waveform version
%-------------------------------------------------------

function states = sync_state_machine2(states, rx_uw)
  ofdm_load_const;
  next_state = states.sync_state;
  states.sync_start = states.sync_end = 0;
  
  if strcmp(states.sync_state,'search') 
    if states.timing_valid
      states.sync_start = 1; states.sync_counter = 0;
      next_state = 'trial';
    end
  end

  states.uw_errors = sum(xor(tx_uw,rx_uw));
 
  if strcmp(states.sync_state,'trial')
    if strcmp(states.sync_state,'trial')
      if states.uw_errors < states.bad_uw_errors;
        next_state = "synced";
        states.frame_count = Nuwframes;
        states.modem_frame = Nuwframes;
      else
        states.sync_counter++;
        if states.sync_counter > Np
          next_state = "search";
        end
      end
    end
  end

  % Note we don't every lose sync, we assume there are a known number of frames being sent,
  % or the packets contain an "end of stream" information.
  if strcmp(states.sync_state,'synced')    
    states.frame_count++;
    states.modem_frame++;
    if (states.modem_frame >= states.Np) states.modem_frame = 0; end
  end
  
  states.last_sync_state = states.sync_state;
  states.sync_state = next_state;
endfunction


% ------------------------------------------------------------------------------
% codec_to_frame_packing - Set up a bunch of constants to support modem frame
%                          construction from LDPC codewords and codec source bits
% ------------------------------------------------------------------------------

function [code_param Nbitspercodecframe Ncodecframespermodemframe] = codec_to_frame_packing(states, mode)
  ofdm_load_const;
  mod_order = 4; bps = 2; modulation = 'QPSK'; mapping = 'gray';

  init_cml('~/cml/');
  if strcmp(mode, "700D")
    load HRA_112_112.txt
    code_param = ldpc_init_user(HRA_112_112, modulation, mod_order, mapping);
    assert(Nbitsperframe == (code_param.coded_bits_per_frame + Nuwbits + Ntxtbits));
    % unused for this mode
    Nbitspercodecframe = Ncodecframespermodemframe = 0;
  end
  if strcmp(mode, "2020")
    load HRA_504_396.txt
    code_param = ldpc_init_user(HRA_504_396, modulation, mod_order, mapping);
    code_param.data_bits_per_frame = 312;
    code_param.coded_bits_per_frame = code_param.data_bits_per_frame + code_param.ldpc_parity_bits_per_frame;
    code_param.coded_syms_per_frame = code_param.coded_bits_per_frame/code_param.bits_per_symbol;
    printf("2020 mode\n");
    printf("ldpc_data_bits_per_frame = %d\n", code_param.ldpc_data_bits_per_frame);
    printf("ldpc_coded_bits_per_frame  = %d\n", code_param.ldpc_coded_bits_per_frame);
    printf("ldpc_parity_bits_per_frame  = %d\n", code_param.ldpc_parity_bits_per_frame);
    printf("data_bits_per_frame = %d\n", code_param.data_bits_per_frame);
    printf("coded_bits_per_frame  = %d\n", code_param.coded_bits_per_frame);
    printf("coded_syms_per_frame  = %d\n", code_param.coded_syms_per_frame);
    printf("ofdm_bits_per_frame  = %d\n", Nbitsperframe);
    Nbitspercodecframe = 52; Ncodecframespermodemframe = 6;
    printf("  Nuwbits: %d  Ntxtbits: %d\n", Nuwbits, Ntxtbits);
    Nparity = code_param.ldpc_parity_bits_per_frame;
    totalbitsperframe = code_param.data_bits_per_frame + Nparity + Nuwbits + Ntxtbits;
    printf("Total bits per frame: %d\n", totalbitsperframe);
    assert(totalbitsperframe == Nbitsperframe);
  end
  if strcmp(mode, "datac1") || strcmp(mode, "datac2") || strcmp(mode, "qam16")
    load H2064_516_sparse.mat
    code_param = ldpc_init_user(HRA, modulation, mod_order, mapping);
  end
  if strcmp(mode, "datac3")
    load H_256_768_22.txt
    code_param = ldpc_init_user(H_256_768_22, modulation, mod_order, mapping);
    Nbitspercodecframe = Ncodecframespermodemframe = -1;
  end
  if strcmp(mode, "datac1") || strcmp(mode, "datac2") || strcmp(mode, "datac3") || strcmp(mode, "qam16")
    printf("ldpc_data_bits_per_frame = %d\n", code_param.ldpc_data_bits_per_frame);
    printf("ldpc_coded_bits_per_frame  = %d\n", code_param.ldpc_coded_bits_per_frame);
    printf("ldpc_parity_bits_per_frame  = %d\n", code_param.ldpc_parity_bits_per_frame);
    printf("Nbitsperpacket  = %d\n", Nbitsperpacket);
    Nparity = code_param.ldpc_parity_bits_per_frame;
    totalbitsperframe = code_param.data_bits_per_frame + Nparity + Nuwbits + Ntxtbits;
    printf("totalbitsperframe = %d\n", totalbitsperframe);
    assert(totalbitsperframe == Nbitsperpacket);
    Nbitspercodecframe = Ncodecframespermodemframe = -1;
  end
endfunction


% ------------------------------------------------------------------------------
% fec_encode - Handle FEC encoding
% ------------------------------------------------------------------------------

function [frame_bits bits_per_frame] = fec_encode(states, code_param, mode, payload_bits, ...
                                                      Ncodecframespermodemframe, Nbitspercodecframe)
  ofdm_load_const;
  if strcmp(mode, "700D") || strcmp(mode, "datac1") || strcmp(mode, "datac2") || strcmp(mode, "datac3") || strcmp(mode, "qam16") 
    frame_bits = LdpcEncode(payload_bits, code_param.H_rows, code_param.P_matrix);
  elseif strcmp(mode, "2020")
    Nunused = code_param.ldpc_data_bits_per_frame - code_param.data_bits_per_frame;
    frame_bits = LdpcEncode([payload_bits zeros(1,Nunused)], code_param.H_rows, code_param.P_matrix);
    % remove unused data bits
    frame_bits = [ frame_bits(1:code_param.data_bits_per_frame) frame_bits(code_param.ldpc_data_bits_per_frame+1:end) ];
  else
    assert(0);
  end
  bits_per_frame = length(frame_bits);
    
endfunction


% test function, kind of like a CRC for QPSK symbols, to compare two vectors

function acc = test_acc(v)
  sre = 0; sim = 0;
  for i=1:length(v)
    x = v(i);
    re = round(real(x)); im = round(imag(x));
    sre += re; sim += im;
    %printf("%d %10f %10f %10f %10f\n", i, re, im, sre, sim);
  end
  acc = sre + j*sim;
end


% Save test bits frame to a text file in the form of a C array
% 
% usage:
%   ofdm_lib; test_bits_ofdm_file
%

function test_bits_ofdm_file
  Ts = 0.018; Tcp = 0.002; Rs = 1/Ts; bps = 2; Nc = 17; Ns = 8;
  states = ofdm_init(bps, Rs, Tcp, Ns, Nc);
  [test_bits_ofdm payload_data_bits codeword] = create_ldpc_test_frame(states);
  printf("%d test bits\n", length(test_bits_ofdm));
  
  f=fopen("../src/test_bits_ofdm.h","wt");
  fprintf(f,"/* Generated by test_bits_ofdm_file() Octave function */\n\n");
  fprintf(f,"const int test_bits_ofdm[]={\n");
  for m=1:length(test_bits_ofdm)-1
    fprintf(f,"  %d,\n",test_bits_ofdm(m));
  endfor
  fprintf(f,"  %d\n};\n",test_bits_ofdm(end));

  fprintf(f,"\nconst int payload_data_bits[]={\n");
  for m=1:length(payload_data_bits)-1
    fprintf(f,"  %d,\n",payload_data_bits(m));
  endfor
  fprintf(f,"  %d\n};\n",payload_data_bits(end));

  fprintf(f,"\nconst int test_codeword[]={\n");
  for m=1:length(codeword)-1
    fprintf(f,"  %d,\n",codeword(m));
  endfor
  fprintf(f,"  %d\n};\n",codeword(end));

  fclose(f);

endfunction

% Get rid of nasty unfiltered stuff either side of OFDM signal
% This may need to be tweaked, or better yet made a function of Nc, if Nc changes
%
% usage:
%  ofdm_lib; make_ofdm_bpf(1);

function bpf_coeff = make_ofdm_bpf(write_c_header_file)
  filt_n = 100;
  Fs = 8000;

  bpf_coeff  = fir2(filt_n,[0 900 1000 2000 2100 4000]/(Fs/2),[0.001 0.001 1 1 0.001 0.001]);

  if write_c_header_file
    figure(1)
    clf;
    h = freqz(bpf_coeff,1,Fs/2);
    plot(20*log10(abs(h)))
    grid minor

    % save coeffs to a C header file

    f=fopen("../src/ofdm_bpf_coeff.h","wt");
    fprintf(f,"/* 1000 - 2000 Hz FIR filter coeffs */\n");
    fprintf(f,"/* Generated by make_ofdm_bpf() in ofdm_lib.m */\n");

    fprintf(f,"\n#define OFDM_BPF_N %d\n\n", filt_n);

    fprintf(f,"float ofdm_bpf_coeff[]={\n");
    for r=1:filt_n
      if r < filt_n
        fprintf(f, "  %f,\n",  bpf_coeff(r));
      else
        fprintf(f, "  %f\n};", bpf_coeff(r));
      end
    end
    fclose(f);
  end

endfunction


