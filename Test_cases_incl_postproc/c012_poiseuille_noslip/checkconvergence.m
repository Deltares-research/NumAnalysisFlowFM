close all;

% Check convergence
time             = nc_varget('dflowfmoutput\poiseuillenoslip_his.nc','time');
u                = nc_varget('dflowfmoutput\poiseuillenoslip_his.nc','x_velocity');
uc               = u(:,10);
for i=2:length(time);
    dudt(i-1)    = (uc(i)-uc(i-1))./(time(i)-time(i-1));
end
clear uc;
semilogy(time(2:end),dudt);
grid on;