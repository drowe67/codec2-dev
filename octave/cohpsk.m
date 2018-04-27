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

    Nd               = sim_in.Nd;     % diveristy
    Ns               = sim_in.Ns;     % step size between pilots
    ldpc_code        = sim_in.ldpc_code;
    rate             = sim_in.ldpc_code_rate; 

    sim_in.bps = bps = 2;

    sim_in.Nsymb         = Nsymb            = framesize/bps;
    sim_in.Nsymbrow      = Nsymbrow         = Nsymb/Nc;
    sim_in.Npilotsframe  = Npilotsframe     = 2;
    sim_in.Nsymbrowpilot = Nsymbrowpilot    = Nsymbrow + Npilotsframe;
    
    if verbose == 2
      printf("Each frame contains %d data bits or %d data symbols, transmitted as %d symbols by %d carriers.", framesize, Nsymb, Nsymbrow, Nc);
      printf("  There are %d pilot symbols in each carrier together at the start of each frame, then %d data symbols.", Npilotsframe, Ns); 
      printf("  Including pilots, the frame is %d symbols long by %d carriers.\n\n", Nsymbrowpilot, Nc);
    end

    sim_in.prev_sym_tx = qpsk_mod([0 0])*ones(1,Nc*Nd);
    sim_in.prev_sym_rx = qpsk_mod([0 0])*ones(1,Nc*Nd);

    sim_in.rx_symb_buf  = zeros(3*Nsymbrow, Nc*Nd);
    sim_in.rx_pilot_buf = zeros(3*Npilotsframe,Nc*Nd);
    sim_in.tx_bits_buf  = zeros(1,2*framesize);

    % pilot sequence is used for phase and amplitude estimation, and frame sync

    pilot = 1 - 2*(rand(Npilotsframe,Nc) > 0.5);
    sim_in.pilot = pilot;
    sim_in.tx_pilot_buf = [pilot; pilot; pilot];
    if sim_in.do_write_pilot_file
      write_pilot_file(pilot, Nsymbrowpilot, Ns, Nsymbrow, Npilotsframe, Nc);
    end

    % we use first 2 pilots of next frame to help with frame sync and fine freq

    sim_in.Nct_sym_buf = 2*Nsymbrowpilot + 2;
    sim_in.ct_symb_buf = zeros(sim_in.Nct_sym_buf, Nc*Nd);

    sim_in.ff_phase = 1;

    sim_in.ct_symb_ff_buf = zeros(Nsymbrowpilot + 2, Nc*Nd);

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

% legacy DQPSK mod for comparative testing

function [tx_symb prev_tx_symb] = bits_to_dqpsk_symbols(sim_in, tx_bits, prev_tx_symb)
    Nc         = sim_in.Nc;
    Nsymbrow   = sim_in.Nsymbrow;

    tx_symb = zeros(Nsymbrow,Nc);

    for c=1:Nc
      for r=1:Nsymbrow
        i = (c-1)*Nsymbrow + r;
        tx_symb(r,c) = qpsk_mod(tx_bits(2*(i-1)+1:2*i));  
        tx_symb(r,c) *= prev_tx_symb(c);
        prev_tx_symb(c) = tx_symb(r,c);
      end
    end
              
endfunction


% legacy DQPSK demod for comparative testing

function [rx_symb rx_bits rx_symb_linear prev_rx_symb] = dqpsk_symbols_to_bits(sim_in, rx_symb, prev_rx_symb)
    Nc         = sim_in.Nc;
    Nsymbrow   = sim_in.Nsymbrow;

    tx_symb = zeros(Nsymbrow,Nc);

    for c=1:Nc
      for r=1:Nsymbrow
        tmp = rx_symb(r,c);
        rx_symb(r,c) *= conj(prev_rx_symb(c))/abs(prev_rx_symb(c));
        prev_rx_symb(c) = tmp;
        i = (c-1)*Nsymbrow + r;
        rx_symb_linear(i) = rx_symb(r,c);
        rx_bits((2*(i-1)+1):(2*i)) = qpsk_demod(rx_symb(r,c));
      end
    end 
              
endfunction


function [tx_symb tx_bits] = bits_to_qpsk_symbols(sim_in, tx_bits, code_param)
    ldpc_code     = sim_in.ldpc_code;
    rate          = sim_in.ldpc_code_rate;
    framesize     = sim_in.framesize;
    Nsymbrow      = sim_in.Nsymbrow;
    Nsymbrowpilot = sim_in.Nsymbrowpilot;
    Nc            = sim_in.Nc;
    Npilotsframe  = sim_in.Npilotsframe;
    Ns            = sim_in.Ns;
    modulation    = sim_in.modulation;
    pilot         = sim_in.pilot;
    Nd            = sim_in.Nd;

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
    
    % insert pilots at start of frame
    
    tx_symb = [pilot(1,:); pilot(2,:); tx_symb;];

    % copy to other carriers (diversity)

    tmp = tx_symb;
    for d=1:Nd-1
      tmp = [tmp tx_symb];
    end
    tx_symb = tmp;

    % ensures energy/symbol is normalised with diversity

    tx_symb = tx_symb/sqrt(Nd);
end


% Symbol rate processing for rx side (demodulator) -------------------------------------------------------

function [rx_symb rx_bits rx_symb_linear amp_ phi_ sig_rms noise_rms cohpsk] = qpsk_symbols_to_bits(cohpsk, ct_symb_buf)
    framesize     = cohpsk.framesize;
    Nsymb         = cohpsk.Nsymb;
    Nsymbrow      = cohpsk.Nsymbrow;
    Nsymbrowpilot = cohpsk.Nsymbrowpilot;
    Nc            = cohpsk.Nc;
    Nd            = cohpsk.Nd;
    Npilotsframe  = cohpsk.Npilotsframe;
    pilot         = cohpsk.pilot;
    verbose       = cohpsk.verbose;
    coh_en        = cohpsk.coh_en;

    % Use pilots to get phase and amplitude estimates We assume there
    % are two samples at the start of each frame and two at the end
    % Note: correlation (averging) method was used initially, but was
    % poor at tracking fast phase changes that we experience on fading
    % channels.  Linear regression (fitting a straight line) works
    % better on fading channels, but increases BER slightly for AWGN
    % channels.

    sampling_points = [1 2 cohpsk.Nsymbrow+3 cohpsk.Nsymbrow+4];
    pilot2 = [cohpsk.pilot(1,:); cohpsk.pilot(2,:); cohpsk.pilot(1,:); cohpsk.pilot(2,:);];
    phi_ = zeros(Nsymbrow, Nc*Nd);
    amp_ = zeros(Nsymbrow, Nc*Nd);
    
    for c=1:Nc*Nd
      %corr = pilot2(:,c)' * ct_symb_buf(sampling_points,c);      
      %phi_(:, c) = angle(corr);
     
      y = ct_symb_buf(sampling_points,c) .* pilot2(:,c-Nc*floor((c-1)/Nc));
      [m b] = linreg(sampling_points, y, length(sampling_points));
      yfit = m*[3 4 5 6] + b;
      phi_(:, c) = angle(yfit);
      %for r=1:Nsymbrow
      %  printf("  %f", phi_(r,c));
      %end
      %printf("\n");
      mag  = sum(abs(ct_symb_buf(sampling_points,c)));
      amp_(:, c) = mag/length(sampling_points);
    end

    % now correct phase of data symbols

    rx_symb = zeros(Nsymbrow, Nc);
    rx_symb_linear = zeros(1, Nsymbrow*Nc*Nd);
    rx_bits = zeros(1, framesize);
    for c=1:Nc*Nd
      for r=1:Nsymbrow
        if coh_en == 1
          rx_symb(r,c) = ct_symb_buf(2+r,c)*exp(-j*phi_(r,c));
        else
          rx_symb(r,c) = ct_symb_buf(2+r,c);
        end
        i = (c-1)*Nsymbrow + r;
        rx_symb_linear(i) = rx_symb(r,c);
        %printf("phi_ %d %d %f %f\n", r,c,real(exp(-j*phi_(r,c))), imag(exp(-j*phi_(r,c))));
      end
    end

    % and finally optional diversity combination and make decn on bits

    for c=1:Nc
      for r=1:Nsymbrow
        i = (c-1)*Nsymbrow + r;
        div_symb = rx_symb(r,c);
        for d=1:Nd-1
          div_symb += rx_symb(r,c + Nc*d);
        end
        rx_bits((2*(i-1)+1):(2*i)) = qpsk_demod(div_symb);
      end
    end

    % Estimate noise power from demodulated symbols.  One method is to
    % calculate the distance of each symbol from the average symbol
    % position. However this is complicated by fading, which means the
    % amplitude of the symbols is constantly changing.
    
    % Now the scatter diagram in a fading channel is a X or cross
    % shape.  The noise can be resolved into two components at right
    % angles to each other.  The component along the the "thickness"
    % of the arms is proportional to the noise power and not affected
    % by fading.  We only use points further along the real axis than
    % the mean amplitude so we keep out of the central nosiey blob
        
    sig_rms = mean(abs(rx_symb_linear));
   
    sum_x = 0;
    sum_xx = 0;
    n = 0;
    for i=1:Nsymb*Nd
      s = rx_symb_linear(i);
      if abs(real(s)) > sig_rms
        sum_x  += imag(s);
        sum_xx += imag(s)*imag(s);
        n++;
      end
    end
   
    noise_var = 0;
    if n > 1
      noise_var = (n*sum_xx - sum_x*sum_x)/(n*(n-1));
    end
    noise_rms = sqrt(noise_var);
      
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

function [spread spread_2ms hf_gain] = init_hf_model(Fs, nsam)

    % convert "spreading" samples from 1kHz carrier at Fss to complex
    % baseband, generated by passing a 1kHz sine wave through PathSim
    % with the ccir-poor model, enabling one path at a time.
    
    Fc = 1000; Fss = 8000;
    fspread = fopen("../raw/sine1k_2Hz_spread.raw","rb");
    spread1k = fread(fspread, "int16")/10000;
    fclose(fspread);
    fspread = fopen("../raw/sine1k_2ms_delay_2Hz_spread.raw","rb");
    spread1k_2ms = fread(fspread, "int16")/10000;
    fclose(fspread);

    % down convert to complex baseband
    spreadbb = spread1k.*exp(-j*(2*pi*Fc/Fss)*(1:length(spread1k))');
    spreadbb_2ms = spread1k_2ms.*exp(-j*(2*pi*Fc/Fss)*(1:length(spread1k_2ms))');

    % remove -2000 Hz image
    b = fir1(50, 5/Fss);
    spread = filter(b,1,spreadbb);
    spread_2ms = filter(b,1,spreadbb_2ms);
   
    % discard first 1000 samples as these were near 0, probably as
    % PathSim states were ramping up

    spread    = spread(1000:length(spread));
    spread_2ms = spread_2ms(1000:length(spread_2ms));

    % change output samples so they are at rate Fs reqd by caller
    
    spread = resample(spread, Fs, Fss);
    spread_2ms = resample(spread_2ms, Fs, Fss);

    % Determine "gain" of HF channel model, so we can normalise
    % carrier power during HF channel sim to calibrate SNR.  I imagine
    % different implementations of ccir-poor would do this in
    % different ways, leading to different BER results.  Oh Well!

    hf_gain = 1.0/sqrt(var(spread(1:nsam))+var(spread_2ms(1:nsam)));
endfunction


function write_pilot_file(pilot, Nsymbrowpilot, Ns, Nsymrow, Npilotsframe, Nc);

  filename = sprintf("../src/cohpsk_defs.h", Npilotsframe, Nc);
  f=fopen(filename,"wt");
  fprintf(f,"/* Generated by write_pilot_file() Octave function */\n\n");
  fprintf(f,"#define NSYMROW      %d   /* number of data symbols on each row (i.e. each carrier) */\n", Nsymrow);
  fprintf(f,"#define NS           %d   /* number of data symbols between pilots                   */\n", Ns);
  fprintf(f,"#define NPILOTSFRAME %d   /* number of pilot symbols on each row                     */\n", Npilotsframe);
  fprintf(f,"#define PILOTS_NC    %d   /* number of carriers                                      */\n\n", Nc);
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

  f=fopen("../src/test_bits_coh.h","wt");
  fprintf(f,"/* Generated by test_bits_coh_file() Octave function */\n\n");
  fprintf(f,"const int test_bits_coh[]={\n");
  for m=1:length(test_bits_coh)-1
    fprintf(f,"  %d,\n",test_bits_coh(m));
  endfor
  fprintf(f,"  %d\n};\n",test_bits_coh(length(test_bits_coh)));
  fclose(f);

endfunction


function [ch_symb rx_timing rx_filt rx_baseband afdmdv f_est] = rate_Fs_rx_processing(afdmdv, ch_fdm_frame, f_est, nsymb, nin, freq_track)
    M = afdmdv.M;
    
    rx_baseband = [];
    rx_filt = [];
    rx_timing = [];

    ch_fdm_frame_index = 1;

    for r=1:nsymb
      % shift signal to nominal baseband, this will put Nc/2 carriers either side of 0 Hz

      [rx_fdm_frame_bb afdmdv.fbb_phase_rx] = freq_shift(ch_fdm_frame(ch_fdm_frame_index:ch_fdm_frame_index + nin - 1), -f_est, afdmdv.Fs, afdmdv.fbb_phase_rx);
      ch_fdm_frame_index += nin;

      % downconvert each FDM carrier to Nc separate baseband signals

      [arx_baseband afdmdv] = fdm_downconvert(afdmdv, rx_fdm_frame_bb, nin);
      [arx_filt afdmdv] = rx_filter(afdmdv, arx_baseband, nin);
      [rx_onesym arx_timing env afdmdv] = rx_est_timing(afdmdv, arx_filt, nin);     

      rx_baseband = [rx_baseband arx_baseband];
      rx_filt     = [rx_filt arx_filt];
      rx_timing    = [rx_timing arx_timing];
      
      ch_symb(r,:) = rx_onesym;

      % we only allow a timing shift on one symbol per frame

      if nin != M
        nin = M;
      end

      % freq tracking, see test_ftrack.m for unit test.  Placed in
      % this function as it needs to work on a symbol by symbol basis
      % rather than frame by frame.  This means the control loop
      % operates at a sample rate of Rs = 50Hz for say 1 Hz/s drift.

      if freq_track
        beta = 0.005;
        g = 0.2;

        % combine difference on phase from last symbol over Nc carriers

        mod_strip = 0;
        for c=1:afdmdv.Nc+1
          adiff = rx_onesym(c) .* conj(afdmdv.prev_rx_symb(c));
          afdmdv.prev_rx_symb(c) = rx_onesym(c);

          % 4th power strips QPSK modulation, by multiplying phase by 4
          % Using the abs value of the real coord was found to help 
          % non-linear issues when noise power was large

          amod_strip = adiff.^4;
          amod_strip = abs(real(amod_strip)) + j*imag(amod_strip);
          mod_strip += amod_strip;
        end
        %plot(mod_strip)
        %printf("modstrip: %f %f\n", real(mod_strip), imag(mod_strip));

        % loop filter made up of 1st order IIR plus integrator.  Integerator
        % was found to be reqd 
        
        afdmdv.filt = (1-beta)*afdmdv.filt + beta*angle(mod_strip);
        f_est += g*afdmdv.filt;
        %printf("filt: %f angle: %f\n", afdmdv.filt, angle(mod_strip));

      end
    end
endfunction


function ct_symb_buf = update_ct_symb_buf(ct_symb_buf, ch_symb, Nct_sym_buf, Nsymbrowpilot)

  % update memory in symbol buffer

  for r=1:Nct_sym_buf-Nsymbrowpilot
    ct_symb_buf(r,:) = ct_symb_buf(r+Nsymbrowpilot,:);
  end
  i = 1;
  for r=Nct_sym_buf-Nsymbrowpilot+1:Nct_sym_buf
    ct_symb_buf(r,:) = ch_symb(i,:);
    i++;
  end
endfunction


% returns index of start of frame and fine freq offset

function [next_sync cohpsk] = frame_sync_fine_freq_est(cohpsk, ch_symb, sync, next_sync)
  ct_symb_buf   = cohpsk.ct_symb_buf;
  Nct_sym_buf   = cohpsk.Nct_sym_buf;
  Rs            = cohpsk.Rs;
  Nsymbrowpilot = cohpsk.Nsymbrowpilot;
  Nc            = cohpsk.Nc;
  Nd            = cohpsk.Nd;

  ct_symb_buf = update_ct_symb_buf(ct_symb_buf, ch_symb, Nct_sym_buf, Nsymbrowpilot);
  cohpsk.ct_symb_buf = ct_symb_buf;
 
  % sample pilots at start of this frame and start of next frame 

  sampling_points = [1 2 cohpsk.Nsymbrow+3 cohpsk.Nsymbrow+4];
  pilot2 = [ cohpsk.pilot(1,:); cohpsk.pilot(2,:); cohpsk.pilot(1,:); cohpsk.pilot(2,:);];

  if sync == 0

    % sample correlation over 2D grid of time and fine freq points

    max_corr = 0;
    for f_fine=-20:0.25:20
%    for f_fine=-1:0.25:1
      f_fine_rect = exp(-j*f_fine*2*pi*sampling_points/Rs)'; % note: this could be pre-computed at init or compile time
      for t=0:cohpsk.Nsymbrowpilot-1
        corr = 0; mag = 0;
        for c=1:Nc*Nd
          f_corr_vec = f_fine_rect .* ct_symb_buf(t+sampling_points,c); % note: this could be pre-computed at init or compile time
          acorr = 0.0;
          for p=1:length(sampling_points)
            acorr += pilot2(p,c-Nc*floor((c-1)/Nc)) * f_corr_vec(p);
            mag   += abs(f_corr_vec(p));
          end
          corr += abs(acorr);
        end
        %printf("  f: %f  t: %d corr: %f mag: %f ratio: %f\n", f_fine, t, corr, mag, corr/mag);
        if corr >= max_corr
          max_corr = corr;
          max_mag = mag;
          cohpsk.ct = t;
          cohpsk.f_fine_est = f_fine;
          cohpsk.ff_rect = exp(-j*f_fine*2*pi/Rs);
        end
      end
    end
    
    printf("  [%d]   fine freq f: %f max_ratio: %f ct: %d\n", cohpsk.frame, cohpsk.f_fine_est, abs(max_corr)/max_mag, cohpsk.ct);
    if abs(max_corr/max_mag) > 0.9
      printf("  [%d]   encouraging sync word! ratio: %f\n", cohpsk.frame, abs(max_corr/max_mag));
      cohpsk.sync_timer = 0;
      next_sync = 1;
    else
      next_sync = 0;
      %printf("  [%d] back to coarse freq offset est...\n", cohpsk.frame) ;
    end
    cohpsk.ratio = abs(max_corr/max_mag);
  end
  
  % single point correlation just to see if we are still in sync

  if sync == 1
    corr = 0; mag = 0;
    f_fine_rect = exp(-j*cohpsk.f_fine_est*2*pi*sampling_points/Rs)';
    for c=1:Nc*Nd
      f_corr_vec = f_fine_rect .* ct_symb_buf(cohpsk.ct+sampling_points,c);
      acorr = 0;
      for p=1:length(sampling_points)
        acorr += pilot2(p, c-Nc*floor((c-1)/Nc)) * f_corr_vec(p);
        mag  += abs(f_corr_vec(p));
      end
      corr += abs(acorr);
    end
    cohpsk.ratio = abs(corr)/mag;
    %printf("f_fine_est: %f ratio: %f\n", cohpsk.f_fine_est, cohpsk.ratio);
  end

endfunction


% misc sync state machine code, just wanted it in a function

function [sync cohpsk] = sync_state_machine(cohpsk, sync, next_sync)

  if sync == 1

    % check that sync is still good, fall out of sync on consecutive bad frames */

    if cohpsk.ratio < 0.8
      cohpsk.sync_timer++;
    else
      cohpsk.sync_timer = 0;            
    end
    % printf("  ratio: %f  sync timer: %d\n", cohpsk.ratio, cohpsk.sync_timer);

    if cohpsk.sync_timer == 10
      printf("  [%d] lost sync ....\n", cohpsk.frame);
      next_sync = 0;
    end
  end

  sync = next_sync;
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
    Nd               = sim_in.Nd;
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
    f_off            = sim_in.f_off;
    div_time_shift   = sim_in.div_timeshift;

    [spread spread_2ms hf_gain] = init_hf_model(Rs, Nsymbrowpilot*(Ntrials+2));

    if strcmp(modulation,'dqpsk')
      Nsymbrowpilot = Nsymbrow;
    end

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

        phase_offset_rect = 1;
        w_offset      = 2*pi*f_off/Rs;
        w_offset_rect = exp(j*w_offset);

        ct_symb_buf = zeros(2*Nsymbrowpilot, Nc*Nd);
        prev_tx_symb = prev_rx_symb = ones(1, Nc*Nd);

        % simulation starts here-----------------------------------
 
        for nn = 1:Ntrials+2
                  
            if ldpc_code
              tx_bits = round(rand(1,framesize*rate));                       
            else
              tx_bits = round(rand(1,framesize));                       
            end

            if strcmp(modulation,'qpsk')

              [tx_symb tx_bits] = bits_to_qpsk_symbols(sim_in, tx_bits, code_param);

              % one frame delay on bits for qpsk

              tx_bits_buf(1:framesize) = tx_bits_buf(framesize+1:2*framesize);
              tx_bits_buf(framesize+1:2*framesize) = tx_bits;

            end
            if strcmp(modulation, 'dqpsk')
              [tx_symb prev_tx_symb] = bits_to_dqpsk_symbols(sim_in, tx_bits, prev_tx_symb);
              tx_bits_buf(1:framesize) = tx_bits;
            end

            s_ch = tx_symb;

            % HF channel simulation  ------------------------------------
            
            hf_fading = ones(1,Nsymb);
            if hf_sim

                % separation between carriers.  Note this effectively
                % under samples at Rs, I dont think this matters.
                % Equivalent to doing freq shift at Fs, then
                % decimating to Rs.

                wsep = 2*pi*(1+0.5);  % e.g. 75Hz spacing at Rs=50Hz, alpha=0.5 filters

                hf_model(hf_n, :) = zeros(1,Nc*Nd);
                
                for r=1:Nsymbrowpilot
                  for c=1:Nd*Nc
                    if c > Nc 
                      time_shift = sim_in.div_timeshift;
                    else
                      time_shift = 1;
                    end
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
              for c=1:Nd*Nc
                 s_ch_tx_log = [s_ch_tx_log s_ch(r,c)];
              end
            end

            % AWGN noise and phase/freq offset channel simulation
            % 0.5 factor ensures var(noise) == variance , i.e. splits power between Re & Im

            noise = sqrt(variance*0.5)*(randn(Nsymbrowpilot,Nc*Nd) + j*randn(Nsymbrowpilot,Nc*Nd));
            noise_log = [noise_log noise];

            for r=1:Nsymbrowpilot
              s_ch(r,:) *= phase_offset_rect;
              phase_offset_rect *= w_offset_rect;
            end
            s_ch += noise;

            ct_symb_buf(1:Nsymbrowpilot,:) = ct_symb_buf(Nsymbrowpilot+1:2*Nsymbrowpilot,:);
            ct_symb_buf(Nsymbrowpilot+1:2*Nsymbrowpilot,:) = s_ch;

            if strcmp(modulation,'qpsk')
              [rx_symb rx_bits rx_symb_linear amp_ phi_ EsNo_ sim_in] = qpsk_symbols_to_bits(sim_in, ct_symb_buf(1:Nsymbrowpilot+Npilotsframe,:));                                 
              phi_log = [phi_log; phi_];
              amp_log = [amp_log; amp_];
            end
            if strcmp(modulation,'dqpsk')
              [rx_symb rx_bits rx_symb_linear prev_rx_symb] = dqpsk_symbols_to_bits(sim_in, s_ch, prev_rx_symb);                                 
            end
                        
            % Wait until we have enough frames to do pilot assisted phase estimation

            if nn > 1
              rx_symb_log = [rx_symb_log rx_symb_linear];
              %EsNo__log = [EsNo__log EsNo_];

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

              printf("EsNo (dB): %3.1f Terrs: %d Tbits: %d BER %5.3f QPSK BER theory %5.3f av_tx_pwr: %3.2f",
                     EsNodB, Terrs, Tbits,
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
          x = 1:Nc*Nd;
          EsNodBSurface = 20*log10(abs(hf_model(y,:))) - 10*log10(variance);
          EsNodBSurface(find(EsNodBSurface < -5)) = -5;
          EsNodBSurface(find(EsNodBSurface > 25)) = 25;
          mesh(x,y,EsNodBSurface);
          grid
          axis([1 Nc*Nd 1 Rs*5 -5 25])
          title('HF Channel Es/No');

          if verbose 
            [m n] = size(hf_model);
            av_hf_pwr = sum(sum(abs(hf_model(:,:)).^2))/(m*n);
            printf("average HF power: %3.2f over %d symbols\n", av_hf_pwr, m*n);
          end

       end

       if strcmp(modulation,'qpsk')
          % set up time axis to include gaps for pilots

         [m1 n1] = size(phi_log);
         phi_x = [];
         phi_x_counter = 1;
         p = Ns;
         for r=1:m1
           if p == Ns
             phi_x_counter += Npilotsframe;
             p = 0;
           end
           p++;
           phi_x = [phi_x phi_x_counter++];        
         end

         phi_x -= Nsymbrowpilot; % account for delay in pilot buffer

         figure(5);
         clf
         subplot(211)
         [m n] = size(phi_log);
         plot(phi_x, phi_log(:,2),'r+;Estimated HF channel phase;')
         if hf_sim
           hold on;
           [m n] = size(hf_model);
           plot(angle(hf_model(1:m,2)),'g;HF channel phase;')
           hold off;
         end
         ylabel('Phase (rads)');
         legend('boxoff');
         axis([1 m -1.1*pi 1.1*pi])

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
         axis([1 m 0 3])
       end

       figure(4)
       clf
       stem(Nerrs_log)
       axis([1 length(Nerrs_log) 0 max(Nerrs_log)+1])
   end

endfunction



function sim_in = standard_init
  sim_in.verbose          = 1;
  sim_in.do_write_pilot_file = 0;
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

  sim_in.Nd            = 1;
endfunction
