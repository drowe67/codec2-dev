% papr_test.m
%
#{

  TODO:
    [ ] measure BER
    [ ] curves with different clipping
    [ ] companding
    [ ] heat map type scatter diagram
#}

M=160;   % length of one symbol
Nc=16;    % number of carriers
N=80000;  % how many samples to simulate
frames = N/M

tx_phase = zeros(frames,Nc);
tx_phase(:, 1:Nc/2) = pi/2*floor((rand(frames,Nc/2)*4));
tx_phase(:, Nc/2+1:Nc) = pi - tx_phase(:,1:Nc/2);

w = 2*pi/M;
tx = [];

% generate OFDM signal
for f=1:frames
  atx = zeros(1,M);
  for c=1:Nc
    atx += exp(j*(0:M-1)*c*w)*exp(j*tx_phase(f,c));
  end
  tx = [tx atx];
end

% channel

rx = tx;
rx(find(rx>10))=10;
% demodulate

rx_symb = zeros(frames,Nc);
for f=1:frames
  st = (f-1)*M+1; en = f*M;
  for c=1:Nc
   rx_symb(f,c) = sum(exp(-j*(0:M-1)*c*w) .* rx(st:en))/M;
  end
end

figure(1); clf;
plot(abs(tx(1:10*M)))
figure(2); clf; hist(abs(rx),25)
figure(3); clf; plot(rx_symb,'+')
printf("PAPR: %f dB\n", 20*log10(max(abs(rx))/mean(abs(rx))))
  
