% bfq19ssa.m
%
% David Rowe Dec 2015
%
% Working for 100mW class A small signal amp using the BFQ19

rfdesign;

w = 2*pi*150E6;

% BFQ19 Vce=10V Ic=50mA 100MHz

S11 = 0.251*exp(j*(-142.7)*pi/180);
S21 = 20.28*exp(j*(103.1)*pi/180);
S12 = 0.030*exp(j*(72.9)*pi/180);
S22 = 0.290*exp(j*(-61.9)*pi/180);

% Lets check stability

Ds = S11*S22-S12*S21;
Knum = 1 + abs(Ds)^2 - abs(S11)^2 - abs(S22)^2;
Kden = 2*abs(S21)*abs(S12);
K = Knum/Kden
figure(1);
clf
scCreate;

if K < 1
  C1 = S11 - Ds*conj(S22);
  C2 = S22 - Ds*conj(S11);
  rs1 = conj(C1)/(abs(S11)^2-abs(Ds)^2);           % centre of input stability circle
  ps1 = abs(S12*S21/(abs(S11)^2-abs(Ds)^2));       % radius of input stability circle
  rs2 = conj(C2)/(abs(S22)^2-abs(Ds)^2);           % centre of input stability circle
  ps2 = abs(S12*S21/(abs(S22)^2-abs(Ds)^2));       % radius of input stability circle

  s(1,1)=S11; s(1,2)=S12; s(2,1)=S21; s(2,2)=S22;
  plotStabilityCircles(s)
end

% determine collector load Rl for our desired power output

P = 0.1;
Irms = 0.02;
Rl = P/(Irms*Irms);

% choose gammaL based on Rl

zo = Rl/50;
[magL,angleL] = ztog(zo);
gammaL = magL*exp(j*angleL*pi/180);

% calculate gammaS and Zi and plot

gammaS = conj(S11 + ((S12*S21*gammaL)/(1 - (gammaL*S22))));
[zi Zi] = gtoz(abs(gammaS), angle(gammaS)*180/pi,50);

scAddPoint(zi);
scAddPoint(zo);

% Design ouput matching network

Ro = 50;
[Xs Xp] = z_match(Ro, Rl)
Cos = 1/(w*Xs);
Lop  = Xp/w;

printf("Output Matching:\n");
printf("  Rl = %3.1f  Ro = %3.1f\n", Rl, Ro);
printf("  Xp = %3.1f Xs = %3.1f\n", Xp, Xs);
printf("  Cos = %3.1f pF Lop = %3.1f nH\n", Cos*1E12, Lop*1E9);

% design input matching network between 50 ohms source and 10 ohms at base

Rb = real(Zi); Rs = 50;

[Xs Xp] = z_match(Rb, Rs);

Lip = Xp/w;
Cis = 1/(w*Xs);

printf("Input Matching:\n");
printf("  Xs = %3.1f Xp = %3.1f\n", Xs, Xp);
printf("  Lp = %3.1f nH Cs = %3.1f pF\n", Lip*1E9, Cis*1E12);
