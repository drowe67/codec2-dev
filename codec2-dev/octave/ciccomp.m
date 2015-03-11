% ciccomp.m
% Brady O'Brien 9 Mar 2015
% CIC Filter compensation helper

graphics_toolkit ("gnuplot");

cicn = 10;                       %delay for CIC filter
N = 10;                         %input interpolation rate
csf = 256;                  %scaling factor for CIC filter conversion
fd = 80e3;                       %DAC frequency
fs = fd/N;                   %Input sampling frequency
fi= fd;                     %freq of DSP
fc1 = fi/4;                     %Frequency of initial upconversion
ciccb=[-0.029626    0.252638   -2.304683   16.332166   -2.304683    0.252638   -0.029626];
%ciccb = ciccb(1:2:length(ciccb))
t = (1:fs);

fdin = zeros(1,length(t));
fdin(1)=1;
fdin=fdin+sin(pi*t*(2000/fs));
%fdin = filter(b,1,fdin);

figure(4)
plot(20*log10(abs(fft(fdin))))
fdcout = zeros(1,length(t));
fdin = int64(fdin*1024);

combout1=0;
combout2=0;
combout3=0;
%combout4=0;
ccbuf1=zeros(1,cicn/N);
ccbuf2=zeros(1,cicn/N);
ccbuf3=zeros(1,cicn/N);
%ccbuf4=zeros(1,cicn);

for ii=1:length(fdin)
	combout1 = fdin(ii) - ccbuf1(end);
	combout2 = combout1 - ccbuf2(end);
	combout3 = combout2 - ccbuf3(end);
	%combout4 = combout3 - ccbuf4(end);
	ccbuf1(2:end)=ccbuf1(1:end-1);
	ccbuf2(2:end)=ccbuf2(1:end-1);
	ccbuf3(2:end)=ccbuf3(1:end-1);
	%ccbuf4(2:end)=ccbuf4(1:end-1);
	ccbuf1(1)=fdin(ii);
	ccbuf2(1)=combout1;
	ccbuf3(1)=combout2;
	%ccbuf4(1)=combout3;
	fdcout(ii)=combout3;
end

intout1=0;
intout2=0;
intout3=0;
%intout4=0;
fdnext = zeros(1,length(t)*N);
fdnext(1:N:length(t)*N) = fdcout; %Interpolate
fdi1 = fdnext;

for ii=1:length(fdnext)
	intout1 = fdnext(ii) + intout1;
	intout2 = intout1 + intout2;
	intout3 = intout2 + intout3;
	%intout4 = intout3 + intout4;
	fdnext(ii)=intout3;
end

figure(5)
plot(20*log10(abs(fft(fdnext))))
