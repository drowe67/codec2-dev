% fsk.m
% David Rowe Nov 2014

% Simulation to test FSK demod
%
% TODO
%   [ ] Code up mod/non-coh demod/AWGN channel simulation
%   [ ] coh demod
%   [ ] Eb/No verses BER curves
%   [ ] test with pre/de-emphahsis impairments
%       + this will introduce delay, use fir filter, group delay
%   [ ] test with hard limiting      
%       + it's FSK, so AM shouldn't matter?  
%       + also make 8-bit fixed point impl easy
%   [ ] channel simulation of HT/FM radio
%       + filtering, varying modulation index
%   [ ] GMSK

rand('state',1); 
randn('state',1);

function sim_out = ber_test(sim_in)
  Fs        = 48000;
  fmark     = 1200;
  fspace    = 2200;
  Rs        = sim_in.Rs;
  Ts        = Fs/Rs;

  framesize = sim_in.framesize;
  EbNodB    = sim_in.EbNodB;

  for ne = 1:length(EbNodB)
    Nerrs = Terrs = Tbits = 0;

    aEbNodB = EbNodB(ne);
    EbNo = 10^(aEbNodB/10);
    variance = Fs/(Rs*EbNo);

    % Modulator -------------------------------

    tx_bits = round(rand(1, framesize));
    tx = zeros(1,framesize*Ts);
    tx_phase = 0;

    for i=1:framesize
      for k=1:Ts
        if tx_bits(i) == 1
          tx_phase += 2*pi*fmark/Fs;
        else
          tx_phase += 2*pi*fspace/Fs;
        end
        tx_phase = tx_phase - floor(tx_phase/(2*pi))*2*pi;
        tx((i-1)*Ts+k) = exp(j*tx_phase);
      end
    end

    % Channel ---------------------------------

    noise = sqrt(variance)*(randn(1,length(tx)));
    rx    = tx + noise;
    %printf("Eb: %f Eb: %f var No: %f EbNo (meas): %f\n", 
    %Eb, var(tx)*Ts/Fs, var(noise)/Fs, (var(tx)*Ts/Fs)/(var(noise)/Fs));
  
    % Demodulator -----------------------------

    % non-coherent fSK demod.  For coherent tracking phase would be
    % hard with a continuous phjase tx as the reference phase would eb
    % a function of all bits to date.

    mark_dc = exp(j*(0:(Ts-1))*2*pi*fmark/Fs);
    space_dc = exp(j*(0:(Ts-1))*2*pi*fspace/Fs);

    rx_bits = zeros(1, framesize);
    for i=1:framesize
      st = (i-1)*Ts+1;
      en = st+Ts-1;
      mark_int(i) = rx(st:en)*mark_dc';
      space_int(i) = rx(st:en)*space_dc';
      rx_bits(i) = abs(mark_int(i)) > abs(space_int(i));
    end
  
    error_positions = xor( rx_bits, tx_bits );
    Nerrs = sum(error_positions);
    Terrs += Nerrs;
    Tbits += length(tx_bits);

    TERvec(ne) = Terrs;
    BERvec(ne) = Terrs/Tbits;

    printf("EbNo (db): %3.2f Terrs: %d BER: %3.2f \n", aEbNodB, Terrs, Terrs/Tbits);
  end

  sim_out.TERvec = TERvec;
  sim_out.BERvec = BERvec;
endfunction

more off;

sim_in.Rs        = 1200;
sim_in.framesize = 1200;
sim_in.EbNodB    = 5:9;

EbNo  = 10 .^ (sim_in.EbNodB/10);
non_coh_theory.BERvec = 0.5*exp(-EbNo/2);
non_coh_sim = ber_test(sim_in);

figure(1); 
clf;
semilogy(sim_in.EbNodB, non_coh_theory.BERvec,'r;FSK non coherent AWGN theory;')
hold on;
semilogy(sim_in.EbNodB, non_coh_sim.BERvec,'g;FSK non coherent AWGN sim;')
hold off;
grid("minor");

