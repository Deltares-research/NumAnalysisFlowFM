function poiseuillepartialslip;

close all;
format long;

% Own settings
nuh           = 0.1;
wallks        = 0.1;
ny            = 20;
B             = 1e3;
L             = 1e4;
deltaH        = -1e-4;

% % Settings 'horvic'
% nuh           = 0.1;
% wallks        = 0.1;
% ny            = 10;
% B             = 100;
% L             = 400;
% ib            = 5e-5;    % verhang gegeven in "unstruc_sysdoc5.ppt", sheet 49
% deltaH        = -ib.*L; 

% Constant settings
kappa         = 0.41;
g             = 9.81;
y0            = wallks/30;
deltay        = B./ny;
dzdx          = -deltaH./L;

% Formulas
alpha         = kappa./log(1+deltay/2/y0);
C0            = g./2.*dzdx;
C1            = nuh./alpha.*sqrt(B.*C0) + C0.*(B./2-deltay./2).^2;
y             = [-B./2:deltay./2:B./2];
uana          = (-C0.*y.^2 + C1)./nuh;

% Read the output data
usim          = nc_varget('dflowfmoutput\poiseuillepartialslip_his.nc','x_velocity');
vsim          = nc_varget('dflowfmoutput\poiseuillepartialslip_his.nc','y_velocity');

% Read the output data
diflon        = (usim(end,ny/2)-uana(ny))/uana(ny)*100;
diflat        = max(abs(vsim(end,:)));
disp(['Difference near centerline is ' ,num2str(diflon,'%3.4f'),'%']);
disp(['Measure for the lateral velocities is ' ,num2str(diflat,'%3.4e')]);

% Plot settings
fs            = 16;
ms            = 20;

% Plot
plot(uana,y); hold on;
plot(usim(end,:),y(2:2:ny*2),'.','markersize',ms); 
hold off;
grid on;
% xlim([0.13  0.17]);
xlim([0.0  0.3]);
ylim([min(y(:)) max(y(:))]);
xlabel(['Streamwise velocity [m/s]'],'fontsize',fs);
ylabel(['Distance to channel centerline [m]'],'fontsize',fs);
% title(['max u = ',num2str(u(ny),'%1.3f'),' m/s   and   min u = ',num2str(u(2),'%1.3f'),' m/s'],'fontsize',fs);
set(gca,'ytick',[-B/2:100:B/2]);
legend('theoretical curve','computed result','location','west');
set(gca,'fontsize',fs);