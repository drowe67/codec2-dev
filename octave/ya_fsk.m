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
endfunction

function tx = yafsk_mod(states,bits)
  Ts = states.Ts;
  Fs = states.Fs;
  fsyms = states.fsyms;
  tx_phase = states.tx_phase;
  %Map bits into symbols
  syms = states.config.txmap(bits);
  tx = zeros(1,Ts*length(syms));
  
  for ii = (1:length(syms))
    cur_sym_f = fsyms(syms(ii));
    for jj = (1:Ts)
        tx_phase = tx_phase + jj*2*pi*cur_sym_f/Fs;
        tx(ii+jj) = exp(j*tx_phase)
    end

    %Correct phase
    if(tx_phase > 2*pi)
        tx_phase = tx_phase - 2*pi;
    elseif(tx_phase < -2*pi);
        tx_phase = tx_phase + 2*pi
    endif
  end
  states.tx_phase = tx_phase;
endfunction

function bits = yafsk_demod(states,stream)
  
endfunction





