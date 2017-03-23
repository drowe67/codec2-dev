% bpsk.m
% David Rowe Mar 2017
%
% Simulate BPSK and DPSK with varoous phase detection schemes

1;

% Differential BPSK detector (m=2)
% ML detection which gains about 0.5dB (m=3, m=4)
%
% Based on JPL publication 89-38 "Multiple Symbol Differential
% Detection of Uncoded and Trellis Coded MPSK" by Divsalar, Simon,
% Shahshahani.  Thanks Johhn Gibbs NN7F for advice.

function rx_symb = dbpsk_demod(m, r)
  tx_set = [1 -1];

  if m == 2

    % regular DBPSK detection

    rx_symb(i) = r(1) * r(2)'/(abs(r(2)));

  else

    % ML DBPSK detection

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
endfunction


function sim_out = run_sim(sim_in)
  Rs = 50;

  Nbits = sim_in.Nbits;
  EbNodB = sim_in.EbNodB;
  dbpsk = sim_in.dbpsk;
  verbose = sim_in.verbose;
  m = sim_in.m;
  phase_est_mem = sim_in.phase_est_mem;
  assert(mod(phase_est_mem,2) == 1, "phase_est_mem must be odd");
  phase_est_delay = floor(phase_est_mem/2);
  printf("phase_est_mem: %d phase_est_delay: %d\n", phase_est_mem, phase_est_delay);

  woffset = 2*pi*sim_in.freq_offset_hz/Rs;

  
  for nn=1:length(EbNodB)
    EbNo = 10 .^ (EbNodB(nn)/10);
    variance = 1/(EbNo/2);
    noise = sqrt(variance)*(0.5*randn(1,Nbits) + j*0.5*randn(1,Nbits));

    tx = zeros(1,Nbits);
    rx = zeros(1,Nbits);
    tx_bits = rand(1,Nbits) > 0.5;
    prev_tx = 1;
    prev_rx = 1;
    r = ones(1,max(4,phase_est_mem));
    phase_amb = 0;
    phase_offset = sim_in.phase_offset

    phase_offset_log = zeros(1,Nbits);
    phase_amb_log = zeros(1,Nbits);
    phase_est_stripped = zeros(1,Nbits);
    phase_est = zeros(1,Nbits);

    Nerrs = 0;
    for i=1:Nbits
      if dbpsk
        tx(i) = prev_tx * exp(j*pi*tx_bits(i));
        prev_tx = tx(i);
      else
        tx(i) = 1 - 2*tx_bits(i);
      end

      rx(i) = tx(i)*exp(j*phase_offset) + noise(i);
      phase_offset_log(i) = phase_offset;
      phase_offset += woffset;
      if (phase_offset > pi)
        phase_offset -= 2*pi;
      end

      r(2:phase_est_mem) = r(1:phase_est_mem-1);
      r(1) = rx(i);

      if dbpsk
        rx_symb(i) = dbpsk_demod(m, r);
      else
        rx_symb(i) = rx(i);
        
        if phase_est_mem 

          if i >= phase_est_mem

            % demod symbol at centre of phase est window

            centre = i - phase_est_delay;

            % modulation strip

            stripped = r(1:phase_est_mem) .^ 2;

            phase_est_stripped(centre) = angle(sum(stripped))/2;

            % determine if phase has jumped from - -> +  
  
            if (phase_est_stripped(centre-1) < -pi/4) && (phase_est_stripped(centre) > pi/4)
              %printf("- -> +\n");
              phase_amb -= pi;
              if (phase_amb < -pi)
                phase_amb += 2*pi;
              end
            end

            % determine if phase has jumped from + -> -    

            if (phase_est_stripped(centre-1) > pi/4) && (phase_est_stripped(centre) < -pi/4)
              %printf("+ -> -\n");
              phase_amb += pi;
              if (phase_amb > pi)
                phase_amb -= 2*pi;
              end
            end

            phase_amb_log(centre) = phase_amb; 
            phase_est(centre) = phase_est_stripped(centre) + phase_amb;

            % keep phase est in range of -pi .. pi to aid plotting

            if (phase_est(centre) < -pi)
              phase_est(centre) += 2*pi;
            end
            if (phase_est(centre) > pi)
              phase_est(centre) -= 2*pi;
            end

            rx_symb(centre) *= exp(-j*phase_est(centre));
          end
        end
      end
    end

    % if using block based phase est, strip off first few and last few
    % symbols where phase ests are invalid
    
    st = phase_est_delay+1;
    %en = Nbits-phase_est_delay;
    en = 50;
    rx_bits = real(rx_symb) < 0;
    errors = xor(tx_bits, rx_bits);
    Nerrs = sum(errors(st:en));                                         
    printf("EbNodB: %3.2f BER: %4.3f Nbits: %d Nerrs: %d\n", EbNodB(nn), Nerrs/(en-st+1), (en-st+1), Nerrs);

    if verbose
      figure(1); clf; 
      plot(rx_symb(st:en),'+');
      axis([-2 2 -2 2]);
      figure(2); clf; 
      plot(phase_offset_log(st:en),'b+;phase offset;', 'markersize', 10);
      hold on;
      plot(phase_amb_log(st:en),'c-;phase amb;');
      plot(phase_est_stripped(st:en),'ro;phase est stripped;');
      plot(phase_est(st:en),'g*;phase est;');
      hold off;
      axis([1 en-st+1 -pi pi])
      legend('boxoff');
      figure(3); clf; 
      plot(errors(st:en),'+');
      axis([1 en-st+1 -0.5 1.5]);
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
  sim_in.Nbits = 1000;
  sim_in.EbNodB = 4;
  sim_in.dbpsk = 0;
  sim_in.m = 2;
  sim_in.phase_est_mem = 11;
  sim_in.verbose = 1;
  sim_in.phase_offset = 0;
  sim_in.freq_offset_hz = 1;

  run_sim(sim_in);
end


format;
more off;
rand('seed',1);
randn('seed',1);

run_single
%run_curves



