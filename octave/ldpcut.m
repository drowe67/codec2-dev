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

  rate = sim_in.rate; 
  framesize = sim_in.framesize;  
  Ntrials = sim_in.Ntrials;
  EsNodBvec = sim_in.EsNodBvec;
  verbose = sim_in.verbose;
  if isfield(sim_in, "hf_en")
    hf_en = sim_in.hf_en;
  else
    hf_en = 0;
  end

  % Start simulation

  mod_order = 4; bps = 2;
  modulation = 'QPSK';
  mapping = 'gray';

  demod_type = 0;
  decoder_type = 0;
  max_iterations = 100;

  code_param = ldpc_init(rate, framesize, modulation, mod_order, mapping);

  % set up HF model

  if hf_en

    % some typical values, or replace with user supplied
    % Note we need to specify a symbol rate Rs

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

    % note we ask for 10% more samples than we need, as
    % doppler_spread() function can sometimes return slightly less
    % than we need due to round off

    spread1 = doppler_spread(dopplerSpreadHz, Rs, Ntrials*framesize*1.1);
    spread2 = doppler_spread(dopplerSpreadHz, Rs, Ntrials*framesize*1.1);
    spread1 = spread1(1:Ntrials*framesize/bps); 
    spread2 = spread2(1:Ntrials*framesize/bps);
    hf_gain = 1.0/sqrt(var(spread1)+var(spread2));
  end

  % ---------------------------------
  % run simulation at each EsNo point
  % ---------------------------------

  for ne = 1:length(EsNodBvec)
    EsNodB = EsNodBvec(ne);
    EsNo = 10^(EsNodB/10);
    variance = 1/EsNo;

    Tbits = Terrs = Ferrs = Terrs_raw = 0;
    
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

    hf_model = ones(1,length(tx_symbols));
    if hf_en

      % Simplified rate Rs simulation model of single carrier.  Note
      % if freq of carrier is 0, model can be simplified further, as
      % path delay term collapses.

      w = 0;      
      hf_model = hf_gain*(spread1 + exp(j*w*path_delay)*spread2);
      rx_symbols = tx_symbols .* abs(hf_model);
    end

    % Add AWGN noise, 0.5 factor ensures var(noise) == variance , i.e. splits power between Re & Im

    noise = sqrt(variance*0.5)*(randn(1,length(tx_symbols)) + j*randn(1,length(tx_symbols)));
    rx_symbols += + noise;

    % Decode a bunch of frames

    for nn = 1: Ntrials        
      st = (nn-1)*code_param.symbols_per_frame + 1;
      en = (nn)*code_param.symbols_per_frame;

      % coded 

      arx_codeword = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, rx_symbols(st:en), EsNo, abs(hf_model(st:en)));
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

      if verbose == 2
        % print "." if frame decoded without errors, 'x' if we can't decode

        if Nerrs > 0, printf('x'),  else printf('.'),  end
        if (rem(nn, 50)==0),  printf('\n'),  end    
      end

      if Nerrs > 0,  Ferrs = Ferrs + 1;  end
      Terrs     += Nerrs;
      Terrs_raw += Nerrs_raw;
      Tbits     += code_param.data_bits_per_frame;        
    end

    if verbose
      printf("\nEsNodB: %3.2f BER: %f Tbits: %d Terrs: %d Ferrs: %d Terrs_raw: %d Raw BER: %3.2f\n", 
             EsNodB, Terrs/Tbits, Tbits, Terrs,  Ferrs, Terrs_raw, Terrs_raw/Tbits)
    end

    sim_out.Tbits(ne) = Tbits;
    sim_out.Terrs(ne) = Terrs;
    sim_out.Ferrs(ne) = Ferrs;
    sim_out.BER(ne)   = Terrs/Tbits;
    sim_out.FER(ne)   = Ferrs/Ntrials;
  end

endfunction


% START SIMULATIONS ---------------------------------------------------------------

more off;
format

% ---------------------------------------------------------------------------------
% 1/ Simplest possible one frame simulation
% ---------------------------------------------------------------------------------

printf("Test 1\n------\n");

mod_order = 4; 
modulation = 'QPSK';
mapping = 'gray';
framesize = 576;         % CML library has a bunch of different framesizes available
rate = 1/2;
demod_type = 0;
decoder_type = 0;
max_iterations = 100;
EsNo = 10;               % decoder needs an estimated channel EsNo (linear ratio, not dB)

code_param = ldpc_init(rate, framesize, modulation, mod_order, mapping);
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

sim_in.rate = 0.5; 
sim_in.framesize = 2*576;  
sim_in.verbose = 2;
sim_in.Ntrials = 100;
sim_in.EsNodBvec = 2;
sim_in.hf_en = 0;
run_simulation(sim_in);


% ---------------------------------------------------------------------------------
% 3/ Lets draw an Eb/No versus BER curve 
% ---------------------------------------------------------------------------------

printf("\n\nTest 3\n------\n");

sim_in.EsNodBvec = -2:10;
sim_in.hf_en = 0;
sim_out = run_simulation(sim_in);
sim_in.hf_en = 1;
sim_out_hf = run_simulation(sim_in);

EbNodB = sim_in.EsNodBvec - 10*log10(2);  % QPSK has two bits/symbol
uncoded_theory = 0.5*erfc(sqrt(10.^(EbNodB/10)));

EbNoLin = 10.^(EbNodB/10);
uncoded_hf_theory = 0.5.*(1-sqrt(EbNoLin./(EbNoLin+1)));

figure(1)
clf
semilogy(EbNodB, uncoded_theory,'r-+;AWGN;')
hold on;
semilogy(EbNodB, uncoded_hf_theory,'r-o;HF;');
semilogy(EbNodB-10*log10(sim_in.rate), sim_out.BER+1E-10,'g-+;AWGN LDPC;');
semilogy(EbNodB-10*log10(sim_in.rate), sim_out_hf.BER+1E-10,'g-o;HF LDPC;');
hold off;
grid('minor')
xlabel('Eb/No (dB)')
ylabel('BER')
axis([min(EbNodB) max(EbNodB) min(uncoded_BER_theory) 1])
legend('boxoff')



