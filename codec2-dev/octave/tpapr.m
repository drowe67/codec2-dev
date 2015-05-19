% tpapr.m
% David Rowe
% 18 May 2015

graphics_toolkit ("gnuplot");
rand('state',1); 

Fs = 8000;
Rs = 50;
Nc = 8;
Fc = 1500;
Fsep = ((1:Nc).^1.2)*75;
%Fsep = (1:Nc)*75 + 5 - 20*rand(1,Nc)
n  = 80000;
t = 1:n;
tx = zeros(1,n);
phi = ones(1,Nc);

figure(1)
clf

for m=1:Nc
    s = cos(phi(m)+t*2*pi*(Rs/2 + Fc + Fsep(m))/Fs);
    tx += s;
    subplot(Nc,1,m);
    plot(s)
end

figure(2)
plot(tx)

tx = tx(length(tx)*0.5:length(tx));
papr = max(tx.*conj(tx)) / mean(tx.*conj(tx));
papr_dB = 10*log10(papr);
printf("PAPR: %4.2f dB\n", papr_dB);
