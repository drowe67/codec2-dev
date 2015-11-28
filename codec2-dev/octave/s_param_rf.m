% s_param_rf.m
%
% David Rowe Nov 2015
%
% RF small signal amplifier design, using equations from "RF Cicruit
% Design" by Chris Bowick

more off;

% BRF92 VCE=5V Ic=5mA 100MHz
 
S11 = 0.727*exp(j*(-43)*pi/180);
S12 = 0.028*exp(j*(69.6)*pi/180);
S21 = 12.49*exp(j*(147)*pi/180);
S22 = 0.891*exp(j*(-16)*pi/180);

% Stability

Ds = S11*S22-S12*S21;
Knum = 1 + abs(Ds)^2 - abs(S11)^2 - abs(S22)^2;
Kden = 2*abs(S21)*abs(S12);
K = Knum/Kden                                    % If > 1 unconditionally stable
                                                 % If < 1 panic
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

% Gain circle

D2 = abs(S22)^2-abs(Ds)^2;
C2 = S22 - Ds*conj(S11);
GdB = 20; Glin = 10^(GdB/10);                     % lets shoot for 20dB gain
G = Glin/(abs(S21)^2);
r0 = G*conj(C2)/(1+D2*G);                                             % centre of gain circle
p0 = sqrt(1 - 2*K*abs(S12*S21)*G + (abs(S12*S21)^2)*(G^2))/(1+D2*G);  % radius of gain circle

scAddCircle(abs(r0),angle(r0)*180/pi,p0,'g')
printf("Green is the %3.1f dB constant gain circle for gammaL\n",GdB);

% Choose a gammaL on the gain circle

gammaL = 0.8 -j*0.4;

% Caclulate gammaS and make sure it's stable by visual inspection
% compared to stability circle.

gammaS = conj(S11 + ((S12*S21*gammaL)/(1 - (gammaL*S22))));
[zo Zo] = gtoz(abs(gammaL), angle(gammaL)*180/pi,50);
[zi Zi] = gtoz(abs(gammaS), angle(gammaS)*180/pi,50);
scAddPoint(zi);
scAddPoint(zo);

% Lets design the z match for the input

  % put input impedance in parallel form

  Zip = zs_to_zp(Zi);

  % first match real part of impedance

  Rs = 50; Rl = real(Zip);
  [Xs Xp] = z_match(Rs,Rl);
  
  % Modify Xp so transistor input sees conjugate match to Zi

  Xp_match = Xp - imag(Zip);

  % Now convert to real component values

  w = 2*pi*150E6;
  Ls = Xs/w;
  Cp = 1/(w*(-Xp_match));

  printf("Input: Zi = %3.1f + %3.1fj ohms\n", real(Zi), imag(Zi));
  printf("       In parallel form Rp = %3.1f Xp = %3.1fj ohms\n", real(Zip), imag(Zip));
  printf("       So for a conjugate match transistor input wants to see:\n         Rp = %3.1f Xp = %3.1fj ohms\n", real(Zip), -imag(Zip));
  printf("       Rs = %3.1f to Rl = %3.1f ohm matching network Xs = %3.1fj Xp = %3.1fj\n", Rs, Rl, Xs, Xp);
  printf("       with conj match to Zi Xs = %3.1fj Xp = %3.1fj\n", Xs, Xp_match);
  printf("       matching components Ls = %5.3f uH Cp = %4.1f pF\n", Ls*1E6, Cp*1E12);

% Now Z match for output

  Lo = -imag(Zo)/w;
  printf("Output: Zo = %3.1f + %3.1fj ohms\n", real(Zo), imag(Zo));
  printf("        So for a conjugate match transistor output wants to see:\n          Rl = %3.1f Xl = %3.1fj ohms\n", real(Zo), -imag(Zo));
  printf("        Which is a series inductor Lo = %5.3f uH\n", Lo*1E6);

% Helper functions -------------------------------------------------

% convert a parallel R/X to a series R/X

function Zs = zp_to_zs(Zp)
  Xp = j*imag(Zp); Rp = real(Zp);
  Zs = Xp*Rp/(Xp+Rp);
endfunction

% convert a series R/X to a parallel R/X

function Zp = zs_to_zp(Zs)
  Xs = imag(Zs); Rs = real(Zs);
  Q = Xs/Rs;       
  Rp = (Q*Q+1)*Rs;
  Xp = Rp/Q;
  Zp = Rp + j*Xp;
endfunction

% Design a Z match network with a parallel and series reactance
% to match between a low and high resistance:
%
%  /--Xs--+---\
%  |      |   |
% Rlow   Xp  Rhigh
%  |      |   |
%  \------+---/
%        

function [Xs Xp] = z_match(Rlow, Rhigh)
  Q = sqrt(Rhigh/Rlow -1);
  Xs = Q*Rlow;
  Xp = -Rhigh/Q; 
endfunction

