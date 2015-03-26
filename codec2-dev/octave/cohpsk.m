% cohpsk.m
% David Rowe Mar 2015
%
% Coherent PSK modem functions, with support for LDPC and DSSS
% (diversity).

1;

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
    if isscalar(symbol) == 0
        printf("only works with scalars\n");
        return;
    end
    bit0 = real(symbol*exp(j*pi/4)) < 0;
    bit1 = imag(symbol*exp(j*pi/4)) < 0;
    two_bits = [bit1 bit0];
endfunction


% init function for symbol rate processing --------------------------------------------------------

function sim_in = symbol_rate_init(sim_in)
    sim_in.Fs = Fs = 8000;

    modulation       = sim_in.modulation;
    verbose          = sim_in.verbose;
    framesize        = sim_in.framesize;
    Ntrials          = sim_in.Ntrials;
    Esvec            = sim_in.Esvec;
    phase_offset     = sim_in.phase_offset;
    w_offset         = sim_in.w_offset;
    plot_scatter     = sim_in.plot_scatter;

    Rs               = sim_in.Rs;
    Nc               = sim_in.Nc;

    hf_sim           = sim_in.hf_sim;
    nhfdelay         = sim_in.hf_delay_ms*Rs/1000;
    hf_mag_only      = sim_in.hf_mag_only;

    Nchip            = sim_in.Nchip;  % spread spectrum factor
    Np               = sim_in.Np;     % number of pilots to use
    Ns               = sim_in.Ns;     % step size between pilots
    ldpc_code        = sim_in.ldpc_code;
    rate             = sim_in.ldpc_code_rate; 

    sim_in.bps = bps = 2;

    sim_in.Nsymb         = Nsymb            = framesize/bps;
    sim_in.Nsymbrow      = Nsymbrow         = Nsymb/Nc;
    sim_in.Npilotsframe  = Npilotsframe     = Nsymbrow/Ns;
    sim_in.Nsymbrowpilot = Nsymbrowpilot    = Nsymbrow + Npilotsframe + 1;

    printf("Each frame is %d bits or %d symbols, transmitted as %d symbols by %d carriers.",
           framesize, Nsymb, Nsymbrow, Nc);
    printf("  There are %d pilot symbols in each carrier, seperated by %d data/parity symbols.",
           Npilotsframe, Ns);
    printf("  Including pilots, the frame is %d symbols long by %d carriers.\n\n", 
           Nsymbrowpilot, Nc);

    assert(Npilotsframe == floor(Nsymbrow/Ns), "Npilotsframe must be an integer");

    sim_in.prev_sym_tx = qpsk_mod([0 0])*ones(1,Nc*Nchip);
    sim_in.prev_sym_rx = qpsk_mod([0 0])*ones(1,Nc*Nchip);

    sim_in.rx_symb_buf  = zeros(3*Nsymbrow, Nc*Nchip);
    sim_in.rx_pilot_buf = zeros(3*Npilotsframe,Nc*Nchip);
    sim_in.tx_bits_buf  = zeros(1,2*framesize);

    % pilot sequence is used for phase and amplitude estimation, and frame sync

    pilot = 1 - 2*(rand(Npilotsframe,Nc) > 0.5);
    sim_in.pilot = pilot;
    sim_in.tx_pilot_buf = [pilot; pilot; pilot];
    if sim_in.do_write_pilot_file
      write_pilot_file(pilot, Nsymbrowpilot, Ns, Np, Nsymbrow, Npilotsframe, Nc);
    end

    % Init LDPC --------------------------------------------------------------------

    if ldpc_code
        % Start CML library

        currentdir = pwd;
        addpath '/home/david/tmp/cml/mat'    % assume the source files stored here
        cd /home/david/tmp/cml
        CmlStartup                           % note that this is not in the cml path!
        cd(currentdir)
  
        % Our LDPC library

        ldpc;

        mod_order = 4; 
        modulation2 = 'QPSK';
        mapping = 'gray';

        sim_in.demod_type = 0;
        sim_in.decoder_type = 0;
        sim_in.max_iterations = 100;

        code_param = ldpc_init(rate, framesize, modulation2, mod_order, mapping);
        code_param.code_bits_per_frame = framesize;
        code_param.symbols_per_frame = framesize/bps;
        sim_in.code_param = code_param;
    else
        sim_in.rate = 1;
        sim_in.code_param = [];
    end
endfunction


% Symbol rate processing for tx side (modulator) -------------------------------------------------------

function [tx_symb tx_bits prev_sym_tx] = bits_to_qpsk_symbols(sim_in, tx_bits, code_param, prev_sym_tx)
    ldpc_code     = sim_in.ldpc_code;
    rate          = sim_in.ldpc_code_rate;
    framesize     = sim_in.framesize;
    Nsymbrow      = sim_in.Nsymbrow;
    Nsymbrowpilot = sim_in.Nsymbrowpilot;
    Nc            = sim_in.Nc;
    Npilotsframe  = sim_in.Npilotsframe;
    Ns            = sim_in.Ns;
    Nchip         = sim_in.Nchip;
    modulation    = sim_in.modulation;
    pilot         = sim_in.pilot;

    if ldpc_code
        [tx_bits, tmp] = ldpc_enc(tx_bits, code_param);
    end

    % modulate --------------------------------------------

    % organise symbols into a Nsymbrow rows by Nc cols
    % data and parity bits are on separate carriers

    tx_symb = zeros(Nsymbrow,Nc);
    
    for c=1:Nc
      for r=1:Nsymbrow
        i = (c-1)*Nsymbrow + r;
        tx_symb(r,c) = qpsk_mod(tx_bits(2*(i-1)+1:2*i));
      end
    end
    %tx_symb = zeros(Nsymbrow,Nc);

    % insert pilots, one every Ns data symbols

    tx_symb_pilot = zeros(Nsymbrowpilot, Nc);
            
    for p=1:Npilotsframe
      tx_symb_pilot((p-1)*(Ns+1)+1,:)          = pilot(p,:);                 % row of pilots
      %printf("%d %d %d %d\n", (p-1)*(Ns+1)+2, p*(Ns+1), (p-1)*Ns+1, p*Ns);
      tx_symb_pilot((p-1)*(Ns+1)+2:p*(Ns+1),:) = tx_symb((p-1)*Ns+1:p*Ns,:); % payload symbols
    end
    tx_symb = tx_symb_pilot;

    % Append extra col of pilots at the start

    tx_symb = [ pilot(1,:);  tx_symb_pilot];

    % Optionally copy to other carriers (spreading)

    for c=Nc+1:Nc:Nc*Nchip
      tx_symb(:,c:c+Nc-1) = tx_symb(:,1:Nc);
    end
            
    % Optionally DQPSK encode
 
    if strcmp(modulation,'dqpsk')
      for c=1:Nc*Nchip
        for r=1:Nsymbrowpilot
          tx_symb(r,c) *= prev_sym_tx(c);
          prev_sym_tx(c) = tx_symb(r,c);
        end
      end               
    end

    % ensures energy/symbol is normalised when spreading

    tx_symb = tx_symb/sqrt(Nchip);
end


% Symbol rate processing for rx side (demodulator) -------------------------------------------------------

function [rx_symb rx_bits rx_symb_linear amp_linear amp_ phi_ EsNo_ prev_sym_rx sim_in] = qpsk_symbols_to_bits(sim_in, s_ch, prev_sym_rx)
    framesize     = sim_in.framesize;
    Nsymb         = sim_in.Nsymb;
    Nsymbrow      = sim_in.Nsymbrow;
    Nsymbrowpilot = sim_in.Nsymbrowpilot;
    Nc            = sim_in.Nc;
    Npilotsframe  = sim_in.Npilotsframe;
    Ns            = sim_in.Ns;
    Np            = sim_in.Np;
    Nchip         = sim_in.Nchip;
    modulation    = sim_in.modulation;
    pilot         = sim_in.pilot;
    rx_symb_buf   = sim_in.rx_symb_buf;
    rx_pilot_buf  = sim_in.rx_pilot_buf;
    tx_pilot_buf  = sim_in.tx_pilot_buf;
    verbose       = sim_in.verbose;

    % demodulate stage 1

    for r=1:Nsymbrowpilot
      for c=1:Nc*Nchip
        rx_symb(r,c) = s_ch(r, c);
        if strcmp(modulation,'dqpsk')
          tmp = rx_symb(r,c);
          rx_symb(r,c) *= conj(prev_sym_rx(c)/abs(prev_sym_rx(c)));
          prev_sym_rx(c) = tmp;
        end
      end
    end
           
    % strip out pilots

    rx_symb_pilot = rx_symb;
    rx_symb = zeros(Nsymbrow, Nc*Nchip);
    rx_pilot = zeros(Npilotsframe, Nc*Nchip);

    for p=1:Npilotsframe
      % printf("%d %d %d %d %d\n", (p-1)*Ns+1, p*Ns, (p-1)*(Ns+1)+2, p*(Ns+1), (p-1)*(Ns+1)+1);
      rx_symb((p-1)*Ns+1:p*Ns,:) = rx_symb_pilot((p-1)*(Ns+1)+2:p*(Ns+1),:);
      rx_pilot(p,:) = rx_symb_pilot((p-1)*(Ns+1)+1,:);
    end

    % buffer three frames of symbols (and pilots) for phase recovery

    rx_symb_buf(1:2*Nsymbrow,:) = rx_symb_buf(Nsymbrow+1:3*Nsymbrow,:);
    rx_symb_buf(2*Nsymbrow+1:3*Nsymbrow,:) = rx_symb;
    rx_pilot_buf(1:2*Npilotsframe,:) = rx_pilot_buf(Npilotsframe+1:3*Npilotsframe,:);
    rx_pilot_buf(2*Npilotsframe+1:3*Npilotsframe,:) = rx_pilot;
    sim_in.rx_symb_buf = rx_symb_buf;
    sim_in.rx_pilot_buf = rx_pilot_buf;

    % pilot assisted phase estimation and correction of middle frame in rx symb buffer

    rx_symb = rx_symb_buf(Nsymbrow+1:2*Nsymbrow,:);
            
    phi_ = zeros(Nsymbrow, Nc*Nchip);
    amp_ = ones(Nsymbrow, Nc*Nchip);

    for c=1:Nc*Nchip

      if verbose > 2
        printf("phi_   : ");
      end

      for r=1:Nsymbrow
        st = Npilotsframe+1+floor((r-1)/Ns) - floor(Np/2) + 1;
        en = st + Np - 1;
        ch_est = tx_pilot_buf(st:en,c)'*rx_pilot_buf(st:en,c)/Np;
        phi_(r,c) = angle(ch_est);
        amp_(r,c) = abs(ch_est);
        %amp_(r,c) = abs(rx_symb(r,c));
        if verbose > 2
          printf("% 4.3f ", phi_(r,c))
        end
        rx_symb(r,c) *= exp(-j*phi_(r,c));
      end

      if verbose > 2
        printf("\nrx_symb: ");
        for r=1:Nsymbrow
          printf("% 4.3f ", angle(rx_symb(r,c)))
        end
        printf("\nindexes: ");
        for r=1:Nsymbrow
          st = Npilotsframe+1+floor((r-1)/Ns) - floor(Np/2) + 1;
          en = st + Np - 1;
          printf("%2d,%2d  ", st,en)
        end
        printf("\npilots : ");
        for p=1:3*Npilotsframe
          printf("% 4.3f ", angle(rx_pilot_buf(p,c)));
        end 
        printf("\n\n");
      end
    end 
    
    % de-spread
            
    for r=1:Nsymbrow
      for c=Nc+1:Nc:Nchip*Nc
        rx_symb(r,1:Nc) = rx_symb(r,1:Nc) + rx_symb(r,c:c+Nc-1);
        amp_(r,1:Nc)    = amp_(r,1:Nc) + amp_(r,c:c+Nc-1);
      end
    end
           
    % demodulate stage 2

    rx_symb_linear = zeros(1,Nsymb);
    amp_linear = zeros(1,Nsymb);
    rx_bits = zeros(1, framesize);
    for c=1:Nc
      for r=1:Nsymbrow
        i = (c-1)*Nsymbrow + r;
        rx_symb_linear(i) = rx_symb(r,c);
        amp_linear(i) = amp_(r,c);
        rx_bits((2*(i-1)+1):(2*i)) = qpsk_demod(rx_symb(r,c));
      end
    end

    % Estimate noise power from demodulated symbols.  One method is to
    % calculate the distance of each symbol from the average symbol
    % position. However this is complicated by fading, which means the
    % amplitude of the symbols is constantly changing.
    
    % Now the scatter diagram in a fading channel is a X shape.  The
    % noise can be resolved into two components at right angles to
    % each other.  The component along the the "thickness" of the arms
    % is proportional to the noise power and not affected by fading.
        
    v = zeros(1,Nsymb);
    for i=1:Nsymb
      s = rx_symb_linear(i);
      if abs(real(s)) > abs(imag(s))
        v(i) = imag(s);
      else
        v(i) = real(s);
      end
      %printf("s: %f %f  v: %f\n", real(s), imag(s), v(i));
    end

    % Note we are only measuring variance in one axis, as other axis is obscured by fading.  We assume
    % that total noise power is shared between both axis so multiply by sqrt(2) to get an estimate of
    % total noise pwr.  Small constant prevents divide by zero errors on start up.

    No_ = var(v)*sqrt(2) + 1E-6;

    % Estimate signal power
    
    Es_ = mean(amp_linear .^ 2);
 
    EsNo_ = Es_/No_;
    %printf("Es_: %f No_: %f  Es/No: %f  Es/No dB: %f\n", Es_, No_, Es_/No_, 10*log10(EsNo_));
  
    % LDPC decoder requires some amplitude normalisation
    % (AGC), was found to break ow.  So we adjust the symbol
    % amplitudes so that they are an averge of 1

    rx_symb_linear /= mean(amp_linear);
    amp_linear /= mean(amp_linear);
    
endfunction


% Compression, John Gibbs pointed out it's best to perform non-linear
% operations on an oversampled signals as they tend to generate
% broadband noise that will be aliased into passband if bandwidth is
% too low

function y = compress(x, power)

  % oversample by a factor of M

  M = 4;
  Ntap = 47;
  n = length(x);  

  b = fir1(Ntap,1/M);
  xM = zeros(1,M*n);
  for i=1:n
    xM(i*M) = M*x(i);
  end
  
  xM = filter(b,1,xM);

  % non linearity

  yM = sign(xM).*(abs(xM) .^ power);

  % decimate by a factor of M

  yM = filter(b,1,yM);
  y  = yM(1:M:n*M);
    
endfunction


% Init HF channel model from stored sample files of spreading signal ----------------------------------

function [spread spread_2ms hf_gain] = init_hf_model(Fs, Rs, nsam)

    % convert "spreading" samples from 1kHz carrier at Fs to complex
    % baseband, generated by passing a 1kHz sine wave through PathSim
    % with the ccir-poor model, enabling one path at a time.
    
    Fc = 1000; M = Fs/Rs;
    fspread = fopen("../raw/sine1k_2Hz_spread.raw","rb");
    spread1k = fread(fspread, "int16")/10000;
    fclose(fspread);
    fspread = fopen("../raw/sine1k_2ms_delay_2Hz_spread.raw","rb");
    spread1k_2ms = fread(fspread, "int16")/10000;
    fclose(fspread);

    % down convert to complex baseband
    spreadbb = spread1k.*exp(-j*(2*pi*Fc/Fs)*(1:length(spread1k))');
    spreadbb_2ms = spread1k_2ms.*exp(-j*(2*pi*Fc/Fs)*(1:length(spread1k_2ms))');

    % remove -2000 Hz image
    b = fir1(50, 5/Fs);
    spread = filter(b,1,spreadbb);
    spread_2ms = filter(b,1,spreadbb_2ms);
   
    % discard first 1000 samples as these were near 0, probably as
    % PathSim states were ramping up

    spread    = spread(1000:length(spread));
    spread_2ms = spread_2ms(1000:length(spread_2ms));

    % decimate down to Rs
    
    spread = spread(1:M:length(spread));
    spread_2ms = spread_2ms(1:M:length(spread_2ms));

    % Determine "gain" of HF channel model, so we can normalise
    % carrier power during HF channel sim to calibrate SNR.  I imagine
    % different implementations of ccir-poor would do this in
    % different ways, leading to different BER results.  Oh Well!

    hf_gain = 1.0/sqrt(var(spread(1:nsam))+var(spread_2ms(1:nsam)));
endfunction


function write_pilot_file(pilot, Nsymbrowpilot, Ns, Np, Nsymrow, Npilotsframe, Nc);

  filename = sprintf("../src/cohpsk_defs.h", Npilotsframe, Nc);
  f=fopen(filename,"wt");
  fprintf(f,"/* Generated by write_pilot_file() Octave function */\n\n");
  fprintf(f,"#define NSYMROW      %d   /* number of data symbols for each row (i.e. each carrier) */\n", Nsymrow);
  fprintf(f,"#define NS           %d   /* number of data symbols between pilots                   */\n", Ns);
  fprintf(f,"#define NP           %d   /* number of pilots to use for channel est                 */\n", Np);
  fprintf(f,"#define NPILOTSFRAME %d   /* number of pilot symbols of each row                     */\n", Npilotsframe);
  fprintf(f,"#define PILOTS_NC    %d   /* number of carriers in coh modem                         */\n\n", Nc);
  fprintf(f,"#define NSYMROWPILOT %d   /* length of row after pilots inserted                     */\n\n", Nsymbrowpilot);
  fclose(f);

  filename = sprintf("../src/pilots_coh.h", Npilotsframe, Nc);
  f=fopen(filename,"wt");
  fprintf(f,"/* Generated by write_pilot_file() Octave function */\n\n");
  fprintf(f,"float pilots_coh[][PILOTS_NC]={\n");
  for r=1:Npilotsframe
    fprintf(f, "  {");
    for c=1:Nc-1
      fprintf(f, "  %f,", pilot(r, c));
    end
    if r < Npilotsframe
      fprintf(f, "  %f},\n", pilot(r, Nc));
    else
      fprintf(f, "  %f}\n};", pilot(r, Nc));
    end
  end
  fclose(f);
endfunction


% Save test bits frame to a text file in the form of a C array

function test_bits_coh_file(test_bits_coh)

  f=fopen("test_bits_coh.h","wt");
  fprintf(f,"/* Generated by test_bits_coh_file() Octave function */\n\n");
  fprintf(f,"const int test_bits_coh[]={\n");
  for m=1:length(test_bits_coh)-1
    fprintf(f,"  %d,\n",test_bits_coh(m));
  endfor
  fprintf(f,"  %d\n};\n",test_bits_coh(length(test_bits_coh)));
  fclose(f);

endfunction


% Frequency offset estimation --------------------------------------------------

function [f_max s_max] = freq_off_est(rx_fdm, tx_pilot, offset, n)

  Fs = 8000;
  nc = 1800;  % portion we wish to correlate over (first 2 rows on pilots)
 
  % downconvert to complex baseband to remove images

  f = 1500;
  foff_rect    = exp(j*2*pi*f*(1:2*n)/Fs);
  tx_pilot_bb  = tx_pilot(1:n) .* foff_rect(1:n)';
  rx_fdm_bb    = rx_fdm(offset:offset+2*n-1) .* foff_rect';

  % remove -2000 Hz image

  b = fir1(50, 1000/Fs);
  tx_pilot_bb_lpf = filter(b,1,tx_pilot_bb);
  rx_fdm_bb_lpf   = filter(b,1,rx_fdm_bb);

  % decimate by M

  M = 4;
  tx_pilot_bb_lpf = tx_pilot_bb_lpf(1:M:length(tx_pilot_bb_lpf));
  rx_fdm_bb_lpf   = rx_fdm_bb_lpf(1:M:length(rx_fdm_bb_lpf));
  n /= M;
  nc /= M;

  % correlate over a range of frequency offsets and delays

  c_max = 0;
  f_n = 1;
  f_range = -75:2.5:75;
  c_log=zeros(n, length(f_range));

  for f=f_range
    foff_rect = exp(j*2*pi*(f*M)*(1:nc)/Fs);
    for s=1:n
      
      c = abs(tx_pilot_bb_lpf(1:nc)' * (rx_fdm_bb_lpf(s:s+nc-1) .* foff_rect'));
      c_log(s,f_n) = c;
      if c > c_max
        c_max = c;
        f_max = f;
        s_max = s;
      end
    end
    f_n++;
    %printf("f: %f c_max: %f f_max: %f s_max: %d\n", f, c_max, f_max, s_max);
  end

  figure(1);
  y = f_range;
  x = max(s_max-25,1):min(s_max+25, n);
  mesh(y,x, c_log(x,:));
  grid
  
  s_max *= M;
  s_max -= floor(s_max/6400)*6400;
  printf("f_max: %f  s_max: %d\n", f_max, s_max);

  % decimated position at sample rate.  need to relate this to symbol
  % rate position.

endfunction


% Rate Rs BER tests ------------------------------------------------------------------------------

function sim_out = ber_test(sim_in)
    sim_in = symbol_rate_init(sim_in);

    Fs               = sim_in.Fs;
    Rs               = sim_in.Rs;
    Ntrials          = sim_in.Ntrials;
    verbose          = sim_in.verbose;
    plot_scatter     = sim_in.plot_scatter;
    framesize        = sim_in.framesize;
    bps              = sim_in.bps;

    Esvec            = sim_in.Esvec;
    ldpc_code        = sim_in.ldpc_code;
    rate             = sim_in.ldpc_code_rate;
    code_param       = sim_in.code_param;
    tx_bits_buf      = sim_in.tx_bits_buf;
    Nsymb            = sim_in.Nsymb;
    Nsymbrow         = sim_in.Nsymbrow;
    Nsymbrowpilot    = sim_in.Nsymbrowpilot;
    Nc               = sim_in.Nc;
    Npilotsframe     = sim_in.Npilotsframe;
    Ns               = sim_in.Ns;
    Np               = sim_in.Np;
    Nchip            = sim_in.Nchip;
    modulation       = sim_in.modulation;
    pilot            = sim_in.pilot;
    prev_sym_tx      = sim_in.prev_sym_tx;
    prev_sym_rx      = sim_in.prev_sym_rx;
    rx_symb_buf      = sim_in.rx_symb_buf;
    tx_pilot_buf     = sim_in.tx_pilot_buf;
    rx_pilot_buf     = sim_in.rx_pilot_buf;

    hf_sim           = sim_in.hf_sim;
    nhfdelay         = sim_in.hf_delay_ms*Rs/1000;
    hf_mag_only      = sim_in.hf_mag_only;

    [spread spread_2ms hf_gain] = init_hf_model(Fs, Rs, Nsymbrowpilot*Ntrials);

    % Start Simulation ----------------------------------------------------------------

    for ne = 1:length(Esvec)
        EsNodB = Esvec(ne);
        EsNo = 10^(EsNodB/10);
    
        variance = 1/EsNo;
        if verbose > 1
            printf("EsNo (dB): %f EsNo: %f variance: %f\n", EsNodB, EsNo, variance);
        end
        
        Terrs = 0;  Tbits = 0;

        s_ch_tx_log      = [];
        rx_symb_log      = [];
        noise_log        = [];
        errors_log       = [];
        Nerrs_log        = [];
        phi_log          = [];
        amp_log          = [];
        EsNo__log        = [];

        ldpc_errors_log = []; ldpc_Nerrs_log = [];

        Terrsldpc = Tbitsldpc = Ferrsldpc = 0;

        % init HF channel

        hf_n = 1;

        phase_offset = 0;
        w_offset     = pi/16;

        % simulation starts here-----------------------------------
 
        for nn = 1:Ntrials+2
                  
            if ldpc_code
              tx_bits = round(rand(1,framesize*rate));                       
            else
              tx_bits = round(rand(1,framesize));                       
            end

            [s_ch tx_bits prev_sym_tx] = bits_to_qpsk_symbols(sim_in, tx_bits, code_param, prev_sym_tx);
   
            tx_bits_buf(1:framesize) = tx_bits_buf(framesize+1:2*framesize);
            tx_bits_buf(framesize+1:2*framesize) = tx_bits;

            % HF channel simulation  ------------------------------------
            
            hf_fading = ones(1,Nsymb);
            if hf_sim

                % separation between carriers.  Note this effectively
                % under samples at Rs, I dont think this matters.
                % Equivalent to doing freq shift at Fs, then
                % decimating to Rs.

                wsep = 2*pi*(1+0.5);  % e.g. 75Hz spacing at Rs=50Hz, alpha=0.5 filters

                hf_model(hf_n, :) = zeros(1,Nc*Nchip);
                
                for r=1:Nsymbrowpilot
                  for c=1:Nchip*Nc
                    time_shift = floor((c-1)*Nsymbrowpilot);
                    ahf_model = hf_gain*(spread(hf_n+time_shift) + exp(-j*c*wsep*nhfdelay)*spread_2ms(hf_n+time_shift));
                    
                    if hf_mag_only
                      s_ch(r,c) *= abs(ahf_model);
                    else
                      s_ch(r,c) *= ahf_model;
                    end
                    hf_model(hf_n, c) = ahf_model;
                  end
                  hf_n++;
                end
            end
           
            % keep a record of each tx symbol so we can check average power

            for r=1:Nsymbrow
              for c=1:Nchip*Nc
                 s_ch_tx_log = [s_ch_tx_log s_ch(r,c)];
              end
            end

            % AWGN noise and phase/freq offset channel simulation
            % 0.5 factor ensures var(noise) == variance , i.e. splits power between Re & Im

            noise = sqrt(variance*0.5)*(randn(Nsymbrowpilot,Nc*Nchip) + j*randn(Nsymbrowpilot,Nc*Nchip));
            noise_log = [noise_log noise];

            s_ch = s_ch + noise;
            
            [rx_symb rx_bits rx_symb_linear amp_linear amp_ phi_ EsNo_ prev_sym_rx sim_in] = qpsk_symbols_to_bits(sim_in, s_ch, prev_sym_rx);                                 

            phi_log = [phi_log; phi_];
            amp_log = [amp_log; amp_];

            % Wait until we have 3 frames to do pilot assisted phase estimation

            if nn > 2 
              rx_symb_log = [rx_symb_log rx_symb_linear];
              EsNo__log = [EsNo__log EsNo_];

              % Measure BER

              error_positions = xor(rx_bits, tx_bits_buf(1:framesize));
              Nerrs = sum(error_positions);
              Terrs += Nerrs;
              Tbits += length(tx_bits);
              errors_log = [errors_log error_positions];
              Nerrs_log = [Nerrs_log Nerrs];

              % Optionally LDPC decode
            
              if ldpc_code
                detected_data = ldpc_dec(code_param, sim_in.max_iterations, sim_in.demod_type, sim_in.decoder_type, ...
                                         rx_symb_linear, min(100,EsNo_), amp_linear);
                error_positions = xor( detected_data(1:framesize*rate), tx_bits_buf(1:framesize*rate) );
                Nerrs = sum(error_positions);
                ldpc_Nerrs_log = [ldpc_Nerrs_log Nerrs];
                ldpc_errors_log = [ldpc_errors_log error_positions];
                if Nerrs
                    Ferrsldpc++;
                end
                Terrsldpc += Nerrs;
                Tbitsldpc += framesize*rate;
              end
            end
          end
           
          TERvec(ne) = Terrs;
          BERvec(ne) = Terrs/Tbits;

            if verbose 
              av_tx_pwr = (s_ch_tx_log * s_ch_tx_log')/length(s_ch_tx_log);

              printf("EsNo est (dB): %3.1f Terrs: %d Tbits: %d BER %4.2f QPSK BER theory %4.2f av_tx_pwr: %3.2f",
                       mean(10*log10(EsNo__log)), Terrs, Tbits,
                       Terrs/Tbits, 0.5*erfc(sqrt(EsNo/2)), av_tx_pwr);
              if ldpc_code
                  printf("\n LDPC: Terrs: %d BER: %4.2f Ferrs: %d FER: %4.2f", 
                         Terrsldpc, Terrsldpc/Tbitsldpc, Ferrsldpc, Ferrsldpc/Ntrials);
              end
              printf("\n");
            end
    end
    
    Ebvec = Esvec - 10*log10(bps);
    sim_out.BERvec          = BERvec;
    sim_out.Ebvec           = Ebvec;
    sim_out.TERvec          = TERvec;
    sim_out.errors_log      = errors_log;
    sim_out.ldpc_errors_log = ldpc_errors_log;

    if plot_scatter
        figure(2);
        clf;
        scat = rx_symb_log .* exp(j*pi/4);
        plot(real(scat), imag(scat),'+');
        title('Scatter plot');
        a = 1.5*max(real(scat)); b = 1.5*max(imag(scat));
        axis([-a a -b b]);

        if hf_sim
          figure(3);
          clf;
        
          y = 1:(hf_n-1);
          x = 1:Nc*Nchip;
          EsNodBSurface = 20*log10(abs(hf_model(y,:))) - 10*log10(variance);
          EsNodBSurface(find(EsNodBSurface < -5)) = -5;
          mesh(x,y,EsNodBSurface);
          grid
          axis([1 (Nc+1)*Nchip 1 Rs*5 -5 15])
          title('HF Channel Es/No');

          if verbose 
            [m n] = size(hf_model);
            av_hf_pwr = sum(sum(abs(hf_model(:,:)).^2))/(m*n);
            printf("average HF power: %3.2f over %d symbols\n", av_hf_pwr, m*n);
          end

       end

        % set up time axis to include gaps for pilots

        [m1 n1] = size(phi_log);
        phi_x = [];
        phi_x_counter = 1;
        p = Ns;
        for r=1:m1
          if p == Ns
            phi_x_counter++;
            p = 0;
          end
          p++;
          phi_x = [phi_x phi_x_counter++];        
        end

        phi_x -= Nsymbrowpilot; % account for delay in pilot buffer

        figure(5);
        clf
        subplot(211)
        plot(phi_x, phi_log(:,2),'r+;Estimated HF channel phase;')
        if hf_sim
          hold on;
          [m n] = size(hf_model);
          plot(angle(hf_model(1:m,2)),'g;HF channel phase;')
          hold off;
        end
        ylabel('Phase (rads)');
        legend('boxoff');

        subplot(212)
        plot(phi_x, amp_log(:,2),'r+;Estimated HF channel amp;')
        if hf_sim
          hold on;
          plot(abs(hf_model(1:m,2)))
          hold off;
        end
        ylabel('Amplitude');
        xlabel('Time (symbols)');
        legend('boxoff');

        figure(4)
        clf
        subplot(211)
        stem(Nerrs_log)
        subplot(212)
        if ldpc_code
          stem(ldpc_Nerrs_log)
        end

   end

endfunction



function sim_in = standard_init
  sim_in.verbose          = 1;
  sim_in.plot_scatter     = 0;

  sim_in.Esvec            = 50; 
  sim_in.Ntrials          = 30;
  sim_in.framesize        = 2;
  sim_in.Rs               = 50;

  sim_in.phase_offset     = 0;
  sim_in.w_offset         = 0;
  sim_in.phase_noise_amp  = 0;

  sim_in.hf_delay_ms      = 2;
  sim_in.hf_sim           = 0;
  sim_in.hf_mag_only      = 0;

  sim_in.Nchip            = 1;
endfunction
