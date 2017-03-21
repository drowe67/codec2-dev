% bpsk.m
% David Rowe Mar 2017
%
% Simulate BPSK and DPSK with varoous phase detection schemes

1;

function sim_out = run_sim(sim_in)
  Rs = 50;

  Nbits = sim_in.Nbits;
  EbNodB = sim_in.EbNodB;
  dbpsk = sim_in.dbpsk;
  verbose = sim_in.verbose;
  m = sim_in.m;
  phase_est_mem = sim_in.phase_est_mem;
  phase_offset = sim_in.phase_offset;
  woffset = 2*pi*sim_in.freq_offset_hz/Rs;

  tx = zeros(1,Nbits);
  rx = zeros(1,Nbits);
  
  for nn=1:length(EbNodB)
    EbNo = 10 .^ (EbNodB(nn)/10);
    variance = 1/(EbNo/2);
    noise = sqrt(variance)*(0.5*randn(1,Nbits) + j*0.5*randn(1,Nbits));

    tx_bits = rand(1,Nbits) > 0.5;
    prev_tx = 1;
    prev_rx = 1;
    r = ones(1,max(4,phase_est_mem));

    Nerrs = 0;
    for i=1:Nbits
      if dbpsk
        tx(i) = prev_tx * exp(j*pi*tx_bits(i));
        prev_tx = tx(i);
      else
        tx(i) = 1 - 2*tx_bits(i);
      end

      rx(i) = tx(i)*exp(j*phase_offset) + noise(i);
      phase_offset += woffset;

      r(2:phase_est_mem) = r(1:phase_est_mem-1);
      r(1) = rx(i);

      if dbpsk
        if m == 2
          %rx_symb(i) = rx(i) * conj(prev_rx)/(abs(prev_rx));
          rx_symb(i) = r(1) * r(2)'/(abs(prev_rx));
          prev_rx = rx(i);
        else
          tx_set = [1 -1];
          max_eta = 0;  rx_symb(i) = tx_set(1);
          for k=1:2
            for k_1=1:2
              for k_2=1:2
                if m == 3
                  eta = abs(r(3) + r(1)*tx_set(k)'*tx_set(k_1)' + r(2)*tx_set(k_1)')^2;
                end
                if m == 4
                  eta = abs(r(4) + r(1)*tx_set(k)'*tx_set(k_1)'*tx_set(k_2)' + r(2)*tx_set(k_1)'*tx_set(k_2)' + r(3)*tx_set(k_2)')^2;
                end
                %printf("  %d %d %f \n", k_1, k, eta);
                if eta > max_eta
                  max_eta = eta;
                  rx_symb(i) = tx_set(k);
                end
              end
            end
          end
        end
      else

        rx_symb(i) = rx(i);
        
        if phase_est_mem 
          % modulation strip
          stripped(i) = angle(sum(r(1:phase_est_mem) .^ 2));
          rx_symb(i) *= exp(-j*stripped(i)/2);
        end
      end
    end

    rx_bits = real(rx_symb) < 0;
    errors = xor(tx_bits, rx_bits);
    Nerrs = sum(errors);
    printf("EbNodB: %3.2f BER: %4.3f Nbits: %d Nerrs: %d\n", EbNodB(nn), Nerrs/Nbits, Nbits, Nerrs);

    if verbose
      figure(1); clf; 
      plot(rx_symb,'+');
      axis([-2 2 -2 2]);
      figure(2); clf; 
      plot(stripped,'+');
      %axis([-2 2 -2 2]);
    end

    sim_out.ber(nn) = Nerrs/Nbits;
  end
endfunction


function run_curves
  sim_in.verbose = 0;
  sim_in.Nbits = 10000;
  sim_in.EbNodB = 0:6;
  sim_in.dbpsk = 0;
  sim_in.m = 2;
  sim_in.phase_est_mem = 0;
  sim_in.phase_offset = 0;

  bpsk_out = run_sim(sim_in);

  sim_in.phase_offset = pi/4;
  sim_in.phase_est_mem = 5;
  bpsk_out_5 = run_sim(sim_in);
  sim_in.phase_est_mem = 10;
  bpsk_out_10 = run_sim(sim_in);
  sim_in.phase_est_mem = 20;
  bpsk_out_20 = run_sim(sim_in);

  figure(3); clf;
  semilogy(sim_in.EbNodB, bpsk_out.ber,'b+-;BPSK;');
  hold on;
  semilogy(sim_in.EbNodB, bpsk_out_5.ber,'g+-;BPSK 5 pt phase est;');
  semilogy(sim_in.EbNodB, bpsk_out_10.ber,'c+-;BPSK 10 pt phase est;');
  semilogy(sim_in.EbNodB, bpsk_out_20.ber,'k+-;BPSK 20 pt phase est;');
  hold off;
  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid;
  legend('boxoff');
  title('Coherent Modn Stripped BPSK');
  print -depsc bpsk_coherent.eps

#{
  sim_in.dbpsk = 1;
  dbpsk_out_2 = run_sim(sim_in);

  sim_in.m = 3;
  dbpsk_out_3 = run_sim(sim_in);
  sim_in.m = 4;
  dbpsk_out_4 = run_sim(sim_in);

  figure(4); clf;
  semilogy(sim_in.EbNodB, bpsk_out.ber,'b+-;BPSK;');
  hold on;
  semilogy(sim_in.EbNodB, dbpsk_out_2.ber,'g+-;DBPSK m=2;');
  semilogy(sim_in.EbNodB, dbpsk_out_3.ber,'c+-;DBPSK m=3;');
  semilogy(sim_in.EbNodB, dbpsk_out_4.ber,'k+-;DBPSK m=4;');
  hold off;
  print -depsc dbpsk_ml.eps

  xlabel('Eb/No (dB)');
  ylabel('BER');
  grid;
  legend('boxoff');
  title('ML DBPSK');
#}
end


function run_single
  sim_in.Nbits = 100;
  sim_in.EbNodB = 40;
  sim_in.dbpsk = 0;
  sim_in.m = 2;
  sim_in.phase_est_mem = 10;
  sim_in.verbose = 1;
  sim_in.phase_offset = pi/4;
  sim_in.freq_offset_hz = .01;

  run_sim(sim_in);
end


format;
more off;
rand('seed',1);
randn('seed',1);

run_single
%run_curves



