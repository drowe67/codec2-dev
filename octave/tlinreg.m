% tlinreg
% David Rowe April 2015
%
% Unit test for linear regression

x = [1 2 7 8];
y = [ -0.70702 + 0.70708i 0.77314 - 0.63442i -0.98083 + 0.19511i 0.99508 - 0.09799i] .* [1 -1 1 -1];
[m b] = linreg(x,y,4);
 
x = [1 2 3 4 5 6 7 8];
yFit = m * x + b;
figure(1)
plot(yFit,'r*')
axis([-1.5 1.5 -1.5 1.5])
