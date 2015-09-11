% fsk4.mf
% 4FSK modem attempt from the DMR spec

graphics_toolkit("gnuplot");

fm;

% Frequency response of the DMR raised cosine filter 
% from ETSI TS 102 361-1 V2.2.1 page 111
fsk4_rcf_resp = @(f) 1.0*(f<=1920) - cos((pi*f)/1920).*1.0.*(f>1920 & f<=2880);

%Maximum positive deviation of amy 4FSK symbol
global fsk4_max_deviation = 1944;

%Deviation of the FSK symbols
global fsk4_symbols = [-1944 -648 648 1944];

function fsk4_states = fsk4_init(fsk4_states,Rs)
    global fsk4_max_deviation;
    global fsk4_symbols;

    Fs = fsk4_states.Fs = 48000;  %Sample rate
    Rs = fsk4_states.Rs = Rs;     %Symbol rate
    M = fsk4_states.M = fsk4_states.Fs/fsk4_states.Rs; %Samples per symbol
    
    %Set up 4FSK raised cosine filter
    rf = (0:(Fs/2));
    tx_filter = fir2(100 ,rf/(Fs/2),fsk4_rcf_resp(rf));
    fsk4_states.tx_filter = tx_filter;
    %Set up the 4FSK symbols
    fsk4_states.symmap = fsk4_symbols / fsk4_max_deviation;
    
    fm_states.Ts = M;
    fm_states.Fs = Fs;
    fm_states.fc = 0;
    fm_states.fm_max = fsk4_max_deviation*2;
    fm_states.fd = fsk4_max_deviation;
    fm_states.pre_emp = fm_states.de_emp = 0;
    fm_states.output_filter = 1;
    fsk4_states.fm_states = analog_fm_init(fm_states);

endfunction 

function [tx, tx_filt] = fsk4_mod(fsk4_states, tx_symbols)
  M = fsk4_states.M;
  nsym = length(tx_symbols);
  nsam = nsym*M;

  tx_stream = zeros(1,nsam);
  for i=1:nsym
    tx_stream(1+(i-1)*M:i*M) = fsk4_states.symmap(tx_symbols(i));
  end
  tx_filt = filter(fsk4_states.tx_filter, 1, tx_stream);
  tx = analog_fm_mod(fsk4_states.fm_states, tx_filt);
endfunction












    
