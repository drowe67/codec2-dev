% yafsk.m
% Yet-another-FSK
% Brady O'Brien 20 October 2015
%
% This is a model for the first attempt at a C FSK modem. Based on fsk_horus and maybe a little fsk4.
% First revision will just be 2fsk
% [x] - Modem init and struct def
% [x] - Direct SDR modulator, probably not FM based
% [.] - Direct SDR non-coherent demodulator
%    [ ] - Core demodulation routine
%    [ ] - Timing offset estimation
%    [ ] - Freq. offset estimation
%    [ ] - Bit slip, maybe
% [ ] - The C port
% [ ] - Some stuff to verify the C port

%clear all;
graphics_toolkit('gnuplot');
fm

%Basic parameters for a simple FSK modem
fsk_setup_info.Rs = 2400;  % Symbol rate
fsk_setup_info.nfsk = 2;      % Number of unique symbols. Must be 2.
fsk_setup_info.Fs = 48000; % Sample frequency
fsk_setup_info.Fsym = fsk_setup_info.Rs; %Symbol spacing
fsk_setup_info.txmap = @(bits) bits+1; %Map TX bits to 2fsk symbols
fsk_setup_info.rxmap = @(syms) syms==2; %Map 2fsk RX symbols to bits


function states = yafsk_init(fsk_config)
  Fs = states.Fs = fsk_config.Fs;
  Rs = states.Rs = fsk_config.Rs;
  nfsk = states.nfsk = fsk_config.nfsk;
  Ts = states.Ts = Fs/Rs;
  Fsym = states.Fsym = fsk_config.Fsym;
  states.config = fsk_config;

  if nfsk != 2
    error("Gotta be 2fsk")
  endif
 
  %Symbol frequencies. Fixed to intervals of Fsym.
  states.fsyms = [-(Fsym/2) (Fsym/2)];
  states.tx_phase = 0;

  states.dc = zeros(1,nfsk);
  states.rx_phi = ones(1,nfsk);
  states.isamp = 0;
  states.sums = zeros(1,nfsk);
endfunction

function [tx states] = yafsk_mod(states,bits)
  Ts = states.Ts;
  Fs = states.Fs;
  fsyms = states.fsyms;
  tx_phase = states.tx_phase;
  %Map bits into symbols
  syms = states.config.txmap(bits);
  tx = zeros(1,Ts*length(syms));
  
  for ii = (1:length(syms))
    cur_sym_f = fsyms(syms(ii));
    tx_phase_i = tx_phase;
    for jj = (1:Ts)
        tx_phase_i = tx_phase + jj*2*pi*cur_sym_f/Fs;
        tx((ii-1)*Ts+jj) = exp(j*tx_phase_i);   
    end
    tx_phase = tx_phase + Ts*2*pi*cur_sym_f/Fs;
    if tx_phase>2*pi
        tx_phase = tx_phase-2*pi;
    elseif tx_phase<-2*pi
        tx_phase = tx_phase+2*pi;
    endif
    %tx_phase_vec = tx_phase + (1:Ts)*2*pi*cur_sym_f/Fs;
    %tx((ii-1)*Ts+1:ii*Ts) = exp(j*tx_phase_vec);
  end
  states.tx_phase = tx_phase;
endfunction

function d = idmp(data, M)
    d = zeros(1,length(data)/M);
    for i = 1:length(d)
      d(i) = sum(data(1+(i-1)*M:i*M));
    end
endfunction

function [bits states] = yafsk_demod_2a(states,rx)
  fine_timing = 1;
  Fs = states.Fs;
  Rs = states.Rs;
  Ts = states.Ts;
  nfsk = states.nfsk;

  phy_f1 = states.rx_phi(1);
  phy_f2 = states.rx_phi(2);

  dphase_f1 = exp(states.fsyms(1)*-j*2*pi/Fs);
  dphase_f2 = exp(states.fsyms(2)*-j*2*pi/Fs);

  sum_f1 = states.sums(1);
  sum_f2 = states.sums(2);

  isamp = states.isamp;
  symcnt = 1;
  syms = [0];
  sums1 = [0];
  sums2 = [0];
  isamp = 1;
  for ii = (1:length(rx))
    phy_f1 *= dphase_f1;   %Spin the oscillators
    phy_f2 *= dphase_f2;

    dcs_f1 = rx(ii)*phy_f1; %Figure out the DC
    dcs_f2 = rx(ii)*phy_f2;

    sum_f1 += dcs_f1; %Integrate
    sum_f2 += dcs_f2;

    isamp += 1;
    if isamp==Ts %If it's time to take a sample and spit out a symbol..
        syms(symcnt) = (abs(sum_f1)>abs(sum_f2))+1; %Spit out a symbol
        symcnt += 1;

        %sums1(symcnt) = abs(sum_f1);
	%sums2(symcnt) = abs(sum_f2);

        sum_f1 = 0;   %Reset integrators
        sum_f2 = 0;
        isamp = 0;    %Reset integrator count
    endif
  end
  plot((1:length(sums1)),sums1,(1:length(sums2)),sums2);

  bits = states.config.rxmap(syms);
endfunction
function [bits states] = yafsk_demod_2(states,rx)
  fine_timing = 1;
  Fs = states.Fs;
  Rs = states.Rs;
  Ts = states.Ts;
  nfsk = states.nfsk;

  rx = rx(fine_timing:length(rx));
  sym_phases = (1:length(rx)).*rot90(states.fsyms)*2*pi/Fs;

  sym_mixers = exp(-j*sym_phases);
  rx_mixd = repmat(rx,nfsk,1).*sym_mixers;
  
  dc1 = abs(idmp(rx_mixd(1,1:length(rx_mixd)),Ts));
  dc2 = abs(idmp(rx_mixd(2,1:length(rx_mixd)),Ts));

  t=(1:length(dc1));
  
  plot(t,dc1,t,dc2)

endfunction





