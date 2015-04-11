
1;

function [m b] = linreg(x,y,n)
  sumx = 0.0;   % sum of x
  sumx2 = 0.0;  % sum of x^2
  sumxy = 0.0;  % sum of x * y
  sumy = 0.0;   % sum of y
  sumy2 = 0.0;  % sum of y**2

  for i=1:n   
    sumx  += x(i);       
    sumx2 += x(i)^2;  
    sumxy += x(i) * y(i);
    sumy  += y(i);      
    sumy2 += y(i)^2; 
  end 

  denom = (n * sumx2 - sumx*sumx);

  if denom == 0
    % singular matrix. can't solve the problem.
    m = 0;
    b = 0;
  else
    m = (n * sumxy  -  sumx * sumy) / denom;
    b = (sumy * sumx2  -  sumx * sumxy) / denom;
  end

endfunction

x = [1 2 7 8];
y = [ -0.70702 + 0.70708i 0.77314 - 0.63442i -0.98083 + 0.19511i 0.99508 - 0.09799i] .* [1 -1 1 -1];
[m b] = linreg(x,y,4);
 
x = [1 2 3 4 5 6 7 8];
yFit = m * x + b;
figure(1)
plot(yFit,'r*')
axis([-1.5 1.5 -1.5 1.5])
