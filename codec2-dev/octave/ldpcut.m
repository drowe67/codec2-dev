% ldpcut.m
%
% David Rowe 18 Dec 2013
%
% Octave LDPC unit test script, based on simulation by Bill Cowley VK5DSP
%

% Start CML library (see CML set up instructions in ldpc.m)

currentdir = pwd;
addpath '/home/david/Desktop/cml/mat' % assume the source files stored here
cd /home/david/Desktop/cml
CmlStartup                            % note that this is not in the cml path!
cd(currentdir)

% Our LDPC library

ldpc;
qpsk;

function sim_out = run_simulation(sim_in)

  % Note this is effective Eb/No of payload data bits, sorta thing we
  % plot on BER versus Eb/No graphs of decoded data.  So if we have a
  % rate 1/2 code, each codeword bit will have Eb/No - 3dB.

  EbNodBvec = sim_in.EbNodBvec;

  % for wimax code frame size specifies code

  if isfield(sim_in, "framesize")
    framesize = sim_in.framesize;
    rate = sim_in.rate; 
  end

  Ntrials = sim_in.Ntrials;
  verbose = sim_in.verbose;
  if isfield(sim_in, "hf_en")
    hf_en = sim_in.hf_en;
  else
    hf_en = 0;
  end
  wimax_en = sim_in.wimax_en;

  % Init LDPC code ------------------------------------

  mod_order = 4; bps = 2;
  modulation = 'QPSK';
  mapping = 'gray';

  demod_type = 0;
  decoder_type = 0;
  max_iterations = 100;

  if wimax_en
    code_param = ldpc_init_wimax(rate, framesize, modulation, mod_order, mapping);
  else 
    load HRA_112_112.txt
    [code_param framesize rate] = ldpc_init_user(HRA_112_112, modulation, mod_order, mapping);
  end

  % set up option HF (multipath) model ------------------------------------

  fading = ones(1,Ntrials*code_param.code_bits_per_frame/bps);

  if hf_en

    % We assume symbols spread acroos Nc OFDM carriers

    Nc = 14; 
    Rs = 50; dopplerSpreadHz = 1.0; path_delay = 1E-3*Rs;

    if isfield(sim_in, "dopplerSpreadHz") 
      dopplerSpreadHz = sim_in.dopplerSpreadHz;
    end
    if isfield(sim_in, "path_delay") 
      path_delay = sim_in.path_delay;
    end
    printf("Doppler Spread: %3.2f Hz Path Delay: %3.2f symbols\n", dopplerSpreadHz, path_delay);

    % reset seed so we get same fading channel on every simulation

    randn('seed',1);

    Ns = Ntrials*code_param.code_bits_per_frame/bps;
    Nr = ceil(Ns/Nc);
    hf_model = zeros(Nr,Nc); 

    % note we ask for 10% more samples than we need, as
    % doppler_spread() function can sometimes return slightly less
    % than we need due to round off

    spread1 = doppler_spread(dopplerSpreadHz, Rs, Nr*1.1);
    spread2 = doppler_spread(dopplerSpreadHz, Rs, Nr*1.1);
    spread1 = spread1(1:Nr); 
    spread2 = spread2(1:Nr);
    hf_gain = 1.0/sqrt(var(spread1)+var(spread2));
  end

  % ----------------------------------
  % run simulation at each EB/No point
  % ----------------------------------

  for ne = 1:length(EbNodBvec)
    randn('seed',1);
    rand('seed',1);

    % Given Eb/No of payload data bits, work out Es/No we need to
    % apply to each channel symbol:
    %
    % i) Each codeword bit gets noise: Eb/No - 3 (for a rate 1/2 code)
    % ii) QPSK means two bits/symbol.: Es/No = Eb/No + 3
    %
    % -> which neatly cancel out ...... (at least for rate 1/2)

    EsNodB = EbNodBvec(ne) + 10*log10(rate) + 10*log10(bps);
    EsNo = 10^(EsNodB/10);
    variance = 1/EsNo;
    hf_r = 1;
    
    Tbits = Terrs = Ferrs = Terrs_raw = Tbits_raw = 0;
    
    tx_bits = [];
    tx_symbols = []; 
    rx_symbols = [];

    % Encode a bunch of frames

    for nn=1:Ntrials        
      atx_bits = round(rand( 1, code_param.data_bits_per_frame));
      tx_bits = [tx_bits atx_bits];
      [tx_codeword atx_symbols] = ldpc_enc(atx_bits, code_param);
      tx_symbols = [tx_symbols atx_symbols];
    end
    
    rx_symbols = tx_symbols;

    % Optional HF (multipath) channel model

    if hf_en

      % Simplified rate Rs OFDM simulation model that doesn't
      % include ISI, just freq filtering.  We assume perfect phase
      % estimation so it's just amplitude distortion.  We distribute
      % symbols across Nc carriers

      Ns = length(tx_symbols);
      w = (1:Nc)*2*pi;  
      rx_symbols = [rx_symbols zeros(1,Nr*Nc-Ns)];  % pad out to integer number of rows

      for r=1:Nr
        for c=1:Nc
          hf_model(hf_r,c) = hf_gain*(spread1(hf_r) + exp(-j*w(c)*path_delay)*spread2(hf_r));
          rx_symbols(Nc*(r-1)+c) *= abs(hf_model(hf_r,c));
          fading(Nc*(r-1)+c) = abs(hf_model(hf_r,c));
        end
        hf_r++;
      end
      rx_symbols = rx_symbols(1:Ns); % remove padding
    end

    % Add AWGN noise, 0.5 factor splits power evenly between Re & Im

    noise = sqrt(variance*0.5)*(randn(1,length(tx_symbols)) + j*randn(1,length(tx_symbols)));
    rx_symbols += noise;

    % Decode a bunch of frames

    for nn = 1: Ntrials        
      st = (nn-1)*code_param.symbols_per_frame + 1;
      en = (nn)*code_param.symbols_per_frame;

      % coded 

      arx_codeword = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, rx_symbols(st:en), EsNo, fading(st:en));
      st = (nn-1)*code_param.data_bits_per_frame + 1;
      en = (nn)*code_param.data_bits_per_frame;
      error_positions = xor(arx_codeword(1:code_param.data_bits_per_frame), tx_bits(st:en));
      Nerrs = sum( error_positions);
        
      % uncoded - to est raw BER compare first half or received frame to tx_bits as code is systematic
      
      rx_bits = [];
      for s=1:code_param.symbols_per_frame*rate
        rx_bits = [rx_bits qpsk_demod(rx_symbols(st+s-1))];
      end
      Nerrs_raw = sum(xor(rx_bits, tx_bits(st:en)));
      Nbits_raw = code_param.data_bits_per_frame;

      if verbose == 2
        % print "." if frame decoded without errors, 'x' if we can't decode

        if Nerrs > 0, printf('x'),  else printf('.'),  end
        if (rem(nn, 50)==0),  printf('\n'),  end    
      end

      if Nerrs > 0,  Ferrs = Ferrs + 1;  end
      Terrs     += Nerrs;
      Tbits     += code_param.data_bits_per_frame;        
      Terrs_raw += Nerrs_raw;
      Tbits_raw += Nbits_raw;
    end

    if verbose
      printf("\nCoded EbNodB: %3.2f BER: %4.3f Tbits: %6d Terrs: %6d FER: %4.3f Tframes: %d Ferrs: %d\n",
             EsNodB, Terrs/Tbits, Tbits, Terrs,  Ferrs/Ntrials, Ntrials, Ferrs);
      EbNodB_raw = EbNodBvec(ne) + 10*log10(rate);
      printf("Raw EbNodB..: %3.2f BER: %4.3f Tbits: %6d Terrs: %6d\n", 
             EbNodB_raw, Terrs_raw/Tbits_raw, Tbits_raw, Terrs_raw);
    end

    sim_out.rate = rate;
    sim_out.Tbits(ne) = Tbits;
    sim_out.Terrs(ne) = Terrs;
    sim_out.Ferrs(ne) = Ferrs;
    sim_out.BER(ne)   = Terrs/Tbits;
    sim_out.FER(ne)   = Ferrs/Ntrials;

    if hf_en && (verbose > 1)
      figure(2); clf;
      plot(real(rx_symbols(Ns/2:Ns)), imag(rx_symbols(Ns/2:Ns)), '+');
      axis([-2 2 -2 2]);
      title('Scatter')

      figure(3); clf;

      % limit mesh plot to Np points to plot quickly
      
      Np = 500;
      step = ceil(hf_r/Np);
      mesh(1:Nc, (1:step:hf_r-1)/Rs, abs(hf_model(1:step:hf_r-1,:)))
      title('HF channel amplitude');
      xlabel('Carrier');
      ylabel('Time (s)');

      figure(4)
      subplot(211); plot(abs(spread1));
      subplot(212); plot(abs(spread2));
      title('HF spreading function amplitudes')
    end
  end

endfunction


% START SIMULATIONS ---------------------------------------------------------------

more off;
format

wimax_en = 1;

% ---------------------------------------------------------------------------------
% 1/ Simplest possible one frame simulation
% ---------------------------------------------------------------------------------

printf("Test 1\n------\n");

mod_order = 4; 
modulation = 'QPSK';
mapping = 'gray';
demod_type = 0;
decoder_type = 0;
max_iterations = 100;

if wimax_en
  framesize = 576*2;         % CML library has a bunch of different framesizes available
  rate = 1/2;
  code_param = ldpc_init_wimax(rate, framesize, modulation, mod_order, mapping);
else
  load HRA_112_112.txt
  [code_param framesize rate] = ldpc_init_user(HRA_112_112, modulation, mod_order, mapping);
end

EsNo = 10;               % decoder needs an estimated channel EsNo (linear ratio, not dB)

tx_bits = round(rand(1, code_param.data_bits_per_frame));
[tx_codeword, qpsk_symbols] = ldpc_enc(tx_bits, code_param);
rx_codeword = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, qpsk_symbols, EsNo);

errors_positions = xor(tx_bits, rx_codeword(1:framesize*rate));
Nerr = sum(errors_positions);
printf("Nerr: %d\n\n", Nerr);


% ---------------------------------------------------------------------------------
% 2/ Run a bunch of trials at just one EsNo point
% ---------------------------------------------------------------------------------

printf("Test 2\n------\n");

sim_in.wimax_en = wimax_en;
if sim_in.wimax_en
  % these are inputs for Wimax mode, e.g. framesize defines code used
  sim_in.rate = 0.5; 
  sim_in.framesize = 576*2;  % long codes smooth over fades but increase latency
end
sim_in.verbose = 2;
sim_in.Ntrials = 1000/4;
sim_in.EbNodBvec = 9;
sim_in.hf_en = 1;
run_simulation(sim_in);


% ---------------------------------------------------------------------------------
% 3/ Lets draw an Eb/No versus BER curve 
% ---------------------------------------------------------------------------------

printf("\n\nTest 3\n------\n");

sim_in.EbNodBvec = -2:10;
sim_in.hf_en = 0;
sim_out = run_simulation(sim_in);
sim_in.hf_en = 1;
sim_out_hf = run_simulation(sim_in);

EbNodB = sim_in.EbNodBvec;
uncoded_theory = 0.5*erfc(sqrt(10.^(EbNodB/10)));

EbNoLin = 10.^(EbNodB/10);
uncoded_hf_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));

figure(1)
clf
semilogy(EbNodB, uncoded_theory,'r-+;AWGN;')
hold on;
semilogy(EbNodB, uncoded_hf_theory,'r-o;HF;');
semilogy(EbNodB, sim_out.BER+1E-10,'g-+;AWGN LDPC;');
semilogy(EbNodB, sim_out_hf.BER+1E-10,'g-o;HF LDPC;');
hold off;
grid('minor')
xlabel('Eb/No (dB)')
ylabel('BER')
axis([min(EbNodB) max(EbNodB) min(uncoded_theory) 1])
legend('boxoff')






