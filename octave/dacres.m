% dacres.m
% David Rowe 18 Feb 2015
% Brady O'Brien 5 Mar 2015
% DAC upconversion simulation

graphics_toolkit ("gnuplot");

cicn = 5;                       %delay for CIC filter
N = 5;                         %input interpolation rate
M = 25;                         %dac interpolation rate
csf = 1024;                  %scaling factor for CIC filter conversion
fd = 2E6;                       %DAC frequency
fs = 2E6/M/N;                   %Input sampling frequency
fi = 2E6/M;                     %freq of first interpolation
fb = 7E5;                       %Bandpass frequency
fc1 = fi/4;                     %Frequency of initial upconversion
ciccb=[-0.029626    0.252638   -2.304683   16.332166   -2.304683    0.252638   -0.029626]; %CIC Compensation FIR
pcicfb = fir1(41,.5); %Interpolation LPF Fir
s1fir = filter(ciccb,1,pcicfb); %Combine compensation and LPF
%s1fir = [-0.00000215, -0.00008715, 0.00073915, -0.00674415, 0.05618415, 0.01629015, -0.19074815, -0.04231615, 0.53620515, 0.09933915, -1.32978715, -0.38797815, 3.97887715, 6.70888315, 3.97887715, -0.38797815, -1.32978715, 0.09933915, 0.53620515, -0.04231615, -0.19074815, ];

t = (1:fs/2);
scpin = e.^(i*(t*pi*2*(3000/8000)));       % initial complex input
%scpin = zeros(1,length(t));
scpin(1) = 1+i;

intstage1 = zeros(1,2*length(scpin));   %First stage of interpolation, 2x
intstage1(1:2:2*length(scpin))=scpin;

scireal = int32(filter(s1fir,1,real(intstage1))*csf);       %separate input into real and imiginary and apply pre-distortion
sciimag = int32(filter(s1fir,1,imag(intstage1))*csf);       % also convert to integer. CIC integrator needs integer to work properly

%Apply 3 stage comb to real
fdin = scireal;
combout1=0;
combout2=0;
combout3=0;
combout4=0;
ccbuf1=zeros(1,cicn/N);
ccbuf2=zeros(1,cicn/N);
ccbuf3=zeros(1,cicn/N);
ccbuf4=zeros(1,cicn/N);

for ii=1:length(fdin)
	combout1 = fdin(ii) - ccbuf1(end);
	combout2 = combout1 - ccbuf2(end);
	combout3 = combout2 - ccbuf3(end);
	combout4 = combout3 - ccbuf4(end);
	ccbuf1(2:end)=ccbuf1(1:end-1);
	ccbuf2(2:end)=ccbuf2(1:end-1);
	ccbuf3(2:end)=ccbuf3(1:end-1);
	ccbuf4(2:end)=ccbuf4(1:end-1);
	ccbuf1(1)=fdin(ii);
	ccbuf2(1)=combout1;
	ccbuf3(1)=combout2;
	ccbuf4(1)=combout3;
	fdcout(ii)=combout4;
end

intout1=0;
intout2=0;
intout3=0;
intout4=0;
fdnext = zeros(1,length(fdcout)*N);
fdnext(1:N:length(fdcout)*N) = fdcout; %Interpolate

for ii=1:length(fdnext)
	intout1 = fdnext(ii) + intout1;
	intout2 = intout1 + intout2;
	intout3 = intout2 + intout3;
	intout4 = intout3 + intout4;
	fdnext(ii)=intout4;
end
scoreal=single(fdnext/(2**20));

fdin=sciimag;

%Apply 3 stage comb to imag
combout1=0;
combout2=0;
combout3=0;
combout4=0;
ccbuf1=zeros(1,cicn/N);
ccbuf2=zeros(1,cicn/N);
ccbuf3=zeros(1,cicn/N);
ccbuf4=zeros(1,cicn/N);

for ii=1:length(fdin)
	combout1 = fdin(ii) - ccbuf1(end);
	combout2 = combout1 - ccbuf2(end);
	combout3 = combout2 - ccbuf3(end);
	combout4 = combout3 - ccbuf4(end);
	ccbuf1(2:end)=ccbuf1(1:end-1);
	ccbuf2(2:end)=ccbuf2(1:end-1);
	ccbuf3(2:end)=ccbuf3(1:end-1);
	ccbuf4(2:end)=ccbuf4(1:end-1);
	ccbuf1(1)=fdin(ii);
	ccbuf2(1)=combout1;
	ccbuf3(1)=combout2;
	ccbuf4(1)=combout3;
	fdcout(ii)=combout4;
end

intout1=0;
intout2=0;
intout3=0;
intout4=0;
fdnext = zeros(1,length(fdcout)*N);
fdnext(1:N:length(fdcout)*N) = fdcout; %Interpolate

for ii=1:length(fdnext)
	intout1 = fdnext(ii) + intout1;
	intout2 = intout1 + intout2;
	intout3 = intout2 + intout3;
	intout4 = intout3 + intout4;
	fdnext(ii)=intout4;
end
scoimag=single(fdnext/(2**20));

%Convert complex to real and shift up to Fs/4
ssout = scoreal.*cos(pi*.5*(1:length(scoreal))) + scoimag.*sin(pi*.5*(1:length(scoimag)));


t = (0:(fi-1));

beta1 = 0.999;
b1x = -2*sqrt(beta1)*cos(2*pi*(fb/fd))
beta2 = beta1 - (1-beta1)*M;

sducin = ssout;
%sducin = cos(pi*.5*t);

sduceq = filter([1 0 beta2],1,sducin);  %pre interpolation notch filter to equalize bandpass after interpolation
sducinterp = zeros(1,length(sduceq)*M);    %interpolate by zero-stuffing
sducinterp(1:M:length(sduceq)*M) = sduceq;
sdac = filter(1,[1 b1x beta1],sducinterp); %select wanted signal

sdac = sdac + median(sdac);  %Center above zero
sdac = sdac / max(sdac);     %normalize
sdac = int32(sdac*2000);     %integerize
sdac = sdac + sdac .^ 5;

figure(1)
subplot(211)
plot(20*log10(abs(fft(sducin)/fi)))
grid
axis([0 fi/2 -20 50])
title('Output from modem');
subplot(212)
plot(20*log10(abs(fft(sduceq/fi))))
grid
axis([0 fi/2 -20 70])
title('After Pre-eq');

figure(2)
subplot(211)
plot(20*log10(abs(fft(sducinterp)/fd)))
grid
axis([0 (fd/2) -20 80])
title('After interpolation');
subplot(212)
plot(20*log10(abs(fft(sdac)/fd)))
grid
title('After bandpass');
%axis([0 (fd/2) -20 80])
