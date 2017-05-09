% ldpcut.m
%
% David Rowe 18 Dec 2013
%
% Octave LDPC unit test script using CML library, based on simulation
% by Bill Cowley VK5DSP

% Libraries we need

ldpc;
qpsk;

function sim_out = run_simulation(sim_in)

  % Note this is effective Eb/No of payload data bits, sorta thing we
  % plot on BER versus Eb/No graphs of decoded data.  So if we have a
  % rate 1/2 code, each codeword bit will have Eb/No - 3dB.

  EbNodBvec = sim_in.EbNodBvec;

  framesize = sim_in.framesize;
  rate = sim_in.rate; 
  Ntrials = sim_in.Ntrials;
  verbose = sim_in.verbose;

  % Init LDPC code ------------------------------------

  mod_order = 4; bps = 2;
  modulation = 'QPSK';
  mapping = 'gray';

  demod_type = 0;
  decoder_type = 0;
  max_iterations = 100;

  code_param = ldpc_init_wimax(rate, framesize, modulation, mod_order, mapping);

  % ----------------------------------
  % run simulation at each Eb/No point
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

    % Add AWGN noise, 0.5 factor splits power evenly between Re & Im

    noise = sqrt(variance*0.5)*(randn(1,length(tx_symbols)) + j*randn(1,length(tx_symbols)));
    rx_symbols += noise;

    % Decode a bunch of frames

    rx_bits_log = [];

    for nn = 1: Ntrials        
      st = (nn-1)*code_param.symbols_per_frame + 1;
      en = (nn)*code_param.symbols_per_frame;

      % coded 

      arx_codeword = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, rx_symbols(st:en), EsNo);
      st = (nn-1)*code_param.data_bits_per_frame + 1;
      en = (nn)*code_param.data_bits_per_frame;
      error_positions = xor(arx_codeword(1:code_param.data_bits_per_frame), tx_bits(st:en));
      Nerrs = sum( error_positions);
      rx_bits_log = [rx_bits_log arx_codeword(1:code_param.data_bits_per_frame)];
        
      % uncoded - to est raw BER compare first half or received frame to tx_bits as code is systematic
      
      raw_rx_bits = [];
      for s=1:code_param.symbols_per_frame*rate
        raw_rx_bits = [raw_rx_bits qpsk_demod(rx_symbols(st+s-1))];
      end
      Nerrs_raw = sum(xor(raw_rx_bits, tx_bits(st:en)));
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
    sim_out.BER(ne)   = Terrs/Tbits;
    sim_out.PER(ne)   = Ferrs/Ntrials;
  end

endfunction

% --------------------------------------------------------------------------------
% START SIMULATIONS
% --------------------------------------------------------------------------------

more off;
format;

% ---------------------------------------------------------------------------------
% Start CML library (see CML set up instructions in ldpc.m)
% ---------------------------------------------------------------------------------

currentdir = pwd;
addpath '/home/david/Desktop/cml/mat' % assume the source files stored here
cd /home/david/Desktop/cml
CmlStartup                            % note that this is not in the cml path!
cd(currentdir)


% ---------------------------------------------------------------------------------
% 1/ Simplest possible one frame simulation
% ---------------------------------------------------------------------------------

printf("\nTest 1\n------\n");

mod_order = 4; 
modulation = 'QPSK';
mapping = 'gray';
demod_type = 0;
decoder_type = 0;
max_iterations = 100;

framesize = 576*2;       % CML library has a bunch of different framesizes available
rate = 1/2;
code_param = ldpc_init_wimax(rate, framesize, modulation, mod_order, mapping);

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

% these are inputs for Wimax mode, e.g. framesize defines code used

sim_in.rate = 0.5; 
sim_in.framesize = 576*4;  % long codes smooth over fades but increase latency

sim_in.verbose = 2;
sim_in.Ntrials = 100;
sim_in.EbNodBvec = 9;
run_simulation(sim_in);


% ---------------------------------------------------------------------------------
% 3/ Lets draw some Eb/No versus BER curves 
% ---------------------------------------------------------------------------------

printf("\n\nTest 3\n------\n");

sim_in.EbNodBvec = -2:10;
sim_out = run_simulation(sim_in);

EbNodB = sim_in.EbNodBvec;
uncoded_awgn_ber_theory = 0.5*erfc(sqrt(10.^(EbNodB/10)));

figure(1)
clf
semilogy(EbNodB, uncoded_awgn_ber_theory,'r-+;AWGN;')
hold on;
semilogy(EbNodB, sim_out.BER+1E-10,'g-+;AWGN LDPC;');
hold off;
grid('minor')
xlabel('Eb/No (dB)')
ylabel('BER')
axis([min(EbNodB) max(EbNodB) 1E-3 1])
legend('boxoff')

