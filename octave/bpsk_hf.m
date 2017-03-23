% bpsk_hf.m
% David Rowe Mar 2017
%
% Rate Rs BPSk simulation to explore phase over multiple carriers in HF channel

1;

function sim_out = run_sim(sim_in)
  Rs = 50;

  Nbits = sim_in.Nbits;
  EbNodB = sim_in.EbNodB;
  verbose = sim_in.verbose;

  for nn=1:length(EbNodB)
    EbNo = 10 .^ (EbNodB(nn)/10);
    variance = 1/(EbNo/2);
    noise = sqrt(variance)*(0.5*randn(1,Nbits) + j*0.5*randn(1,Nbits));

    tx = zeros(1,Nbits);
    rx = zeros(1,Nbits);
    tx_bits = rand(1,Nbits) > 0.5;
    tx = 1 - 2*tx_bits;
    rx = tx + noise;

    rx_bits = real(rx) < 0;
    errors = xor(tx_bits, rx_bits);
    Nerrs = sum(errors);                                         
    printf("EbNodB: %3.2f BER: %4.3f Nbits: %d Nerrs: %d\n", EbNodB(nn), Nerrs/Nbits, Nbits, Nerrs);

    if verbose
      figure(1); clf; 
      plot(rx,'+');
      axis([-2 2 -2 2]);
    end

    sim_out.ber(nn) = Nerrs/Nbits;
  end
endfunction


function run_curves
  sim_in.verbose = 0;
  sim_in.Nbits = 10000;
  sim_in.EbNodB = 0:6;

  bpsk_out = run_sim(sim_in);

  figure(3); clf;
  semilogy(sim_in.EbNodB, bpsk_out.ber,'b+-;BPSK;');
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid;
  legend('boxoff');
  title('Coherent BPSK');

end


function run_single
  sim_in.Nbits = 1000;
  sim_in.EbNodB = 4;
  sim_in.verbose = 1;

  run_sim(sim_in);
end


format;
more off;
rand('seed',1);
randn('seed',1);

run_single
%run_curves



