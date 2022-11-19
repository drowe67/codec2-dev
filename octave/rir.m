function [h]=rir(fs, mic, n, r, rm, src);
%RIR   Room Impulse Response.
%   [h] = RIR(FS, MIC, N, R, RM, SRC) performs a room impulse
%         response calculation by means of the mirror image method.
%
%      FS  = sample rate.
%      MIC = row vector giving the x,y,z coordinates of
%            the microphone.  
%      N   = The program will account for (2*N+1)^3 virtual sources 
%      R   = reflection coefficient for the walls, in general -1<R<1.
%      RM  = row vector giving the dimensions of the room.  
%      SRC = row vector giving the x,y,z coordinates of 
%            the sound source.
%
%   EXAMPLE:
%
%      >>fs=44100;
%      >>mic=[19 18 1.6];
%      >>n=12;
%      >>r=0.3;
%      >>rm=[20 19 21];
%      >>src=[5 2 1];
%      >>h=rir(fs, mic, n, r, rm, src);
%
%   NOTES:
%
%   1) All distances are in meters.
%   2) The output is scaled such that the largest value of the 
%      absolute value of the output vector is equal to one.
%   3) To implement this filter, you will need to do a fast 
%      convolution.  The program FCONV.m will do this. It can be 
%      found on the Mathworks File Exchange at
%      www.mathworks.com/matlabcentral/fileexchange/.  It can also 
%      be found at http://www.sgm-audio.com/research/rir/fconv.m
%   4) A paper has been written on this model.  It is available at:
%      http://www.sgm-audio.com/research/rir/rir.html
%      
%
%Version 3.4.2
%Copyright © 2003 Stephen G. McGovern
%Some of the following comments are references to equations the my paper.
nn=-n:1:n;                            % Index for the sequence
rms=nn+0.5-0.5*(-1).^nn;              % Part of equations 2,3,& 4
srcs=(-1).^(nn);                      % part of equations 2,3,& 4
xi=srcs*src(1)+rms*rm(1)-mic(1);      % Equation 2 
yj=srcs*src(2)+rms*rm(2)-mic(2);      % Equation 3 
zk=srcs*src(3)+rms*rm(3)-mic(3);      % Equation 4 
[i,j,k]=meshgrid(xi,yj,zk);           % convert vectors to 3D matrices
d=sqrt(i.^2+j.^2+k.^2);               % Equation 5
time=round(fs*d/343)+1;               % Similar to Equation 6
              
[e,f,g]=meshgrid(nn, nn, nn);         % convert vectors to 3D matrices
c=r.^(abs(e)+abs(f)+abs(g));          % Equation 9
e=c./d;                               % Equivalent to Equation 10
h=full(sparse(time(:),1,e(:)));       % Equivalent to equation 11
h=h/max(abs(h));                      % Scale output

