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

  % Start simulation

  mod_order = 4; 
  modulation = 'QPSK';
  mapping = 'gray';

  demod_type = 0;
  decoder_type = 0;
  max_iterations = 100;

  code_param = ldpc_init(rate, framesize, modulation, mod_order, mapping);

  for ne = 1:length(EsNodBvec)
    EsNodB = EsNodBvec(ne);
    EsNo = 10^(EsNodB/10);
    variance = 1/EsNo;

    Tbits = Terrs = Ferrs = 0;
    
    tx_bits = [];
    tx_symbols = []; 
    rx_symbols = [];

    % Encode a bunch of frames

    for nn = 1: Ntrials        
      atx_bits = round(rand( 1, code_param.data_bits_per_frame));
      tx_bits = [tx_bits atx_bits];
      [tx_codeword atx_symbols] = ldpc_enc(atx_bits, code_param);
      tx_symbols = [tx_symbols atx_symbols];
    end

    % Add AWGN noise, 0.5 factor ensures var(noise) == variance , i.e. splits power between Re & Im

    noise = sqrt(variance*0.5)*(randn(1,length(tx_symbols)) + j*randn(1,length(tx_symbols)));
    rx_symbols = tx_symbols + noise;

    % Decode a bunch of frames

    for nn = 1: Ntrials        
      st = (nn-1)*code_param.symbols_per_frame + 1;
      en = (nn)*code_param.symbols_per_frame;
      arx_codeword = ldpc_dec(code_param, max_iterations, demod_type, decoder_type, rx_symbols(st:en), EsNo, ones(1,code_param.symbols_per_frame));
      st = (nn-1)*code_param.data_bits_per_frame + 1;
      en = (nn)*code_param.data_bits_per_frame;
      error_positions = xor(arx_codeword(1:code_param.data_bits_per_frame), tx_bits(st:en));
      Nerrs = sum( error_positions);
        
      if verbose == 2
        % print "." if frame decoded without errors, 'x' if we can't decode

        if Nerrs > 0, printf('x'),  else printf('.'),  end
        if (rem(nn, 50)==0),  printf('\n'),  end    
      end

      if Nerrs > 0,  Ferrs = Ferrs + 1;  end
      Terrs = Terrs + Nerrs;
      Tbits = Tbits + code_param.data_bits_per_frame;        
    end

    if verbose
      printf("\nEsNodB: %3.2f BER: %f Tbits: %d Terrs: %d Ferrs: %d\n", EsNodB, Terrs/Tbits, Tbits, Terrs,  Ferrs)
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

% 1/ Simplest possible one frame simulation ---------------------------------------

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

% 2/ Run a bunch of trials at just one EsNo point --------------------------------------

printf("Test 2\n------\n");

sim_in.rate = 0.5; 
sim_in.framesize = 2*576;  
sim_in.verbose = 2;
sim_in.Ntrials = 100;
sim_in.EsNodBvec = 5;
run_simulation(sim_in);

% 3/ Lets draw a Eb/No versus BER curve -------------------------------------------------

printf("\n\nTest 3\n------\n");

sim_in.EsNodBvec = -2:10;
sim_out = run_simulation(sim_in);

EbNodB = sim_in.EsNodBvec - 10*log10(2);  % QPSK has two bits/symbol
uncoded_BER_theory = 0.5*erfc(sqrt(10.^(EbNodB/10)));

figure(1)
clf
semilogy(EbNodB, uncoded_BER_theory,'r;uncoded QPSK theory;')
hold on;
semilogy(EbNodB-10*log10(sim_in.rate), sim_out.BER+1E-10,'g;LDPC coded QPSK simulation;');
hold off;
grid('minor')
xlabel('Eb/No (dB)')
ylabel('BER')
axis([min(EbNodB) max(EbNodB) min(uncoded_BER_theory) 1])
