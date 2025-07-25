function initialSetup;

fclose all; clear all; close all; clc;

% Read file
xk             = nc_varget('planar1d_net.nc','NetNode_x');
yk             = nc_varget('planar1d_net.nc','NetNode_y');

% Define bottom
L              =  300.0;
r0             =  120.0;
dep            =  10;
zk             = -dep.*(1 - (xk-L/2).^2./r0.^2);

% Second input
eta            = 0.23;
psi            = eta.*r0;
samp           = eta.*dep./r0;
hk             = samp.*(2.*(xk-L/2) - psi);

% Put zeros
j              = find(hk<=zk | hk==0);
hk(j)          = -999;

% Assign variables
iniwatk        = [xk yk hk];  
inibott        = [xk yk zk]; 
iniwatk        = sortrows(iniwatk);
inibott        = sortrows(inibott);

% Write files
dlmwrite('iniwatk.xyz', iniwatk, 'delimiter', '\t','precision','%7.7f');
dlmwrite('inibott.xyb', inibott, 'delimiter', '\t','precision','%7.7f');