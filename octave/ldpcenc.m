%test batch file   for the FSO-OOK simulations
% 
% this version uses WiMax eIRA codes and includes erasures 

currentdir = pwd;
thiscomp = computer;

addpath '/home/david/tmp/cml/mat'    % assume the source files stored here
cd /home/david/tmp/cml
CmlStartup              % note that this is not in the cml path!
disp('added cluster path and run CmlStartup')

cd(currentdir)

ldpc;

%
%sim_in.Eprob = 0.1; 
%disp([' test with erasure probability of ' num2str(sim_in.Eprob)]); 

%sim_in.comment = 'test ldpc';
%sim_in.Esvec = 8:1/2:11; 

%sim_in.rate = 3/4;
%sim_in.framesize = 16200

% sim_in.rate = 0.5;
% sim_in.framesize = 204;

sim_in.rate = 3/4; 
sim_in.framesize = 576;  

sim_in.mod_order = 4; 
sim_in.modulation = 'QPSK';
sim_in.mapping = 'gray';

sim_in.Lim_Ferrs= 30;
sim_in.Ntrials =  1000;


if exist('deb')~=1, deb = 0; end
sim_in.deb =    deb;

% init enc and dec
% generate bits
% encode
% decode
% measure BER
ldpc_proc(sim_in, 'test.mat');

