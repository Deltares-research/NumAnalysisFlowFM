function showgridconvergence;

fclose all;
close all;
clc;

% Print?
printen     = 1;

% Check steady state
t           = nc_varget('dflowfmoutput/flatbed_quadrilateral_his.nc','time');
h           = nc_varget('dflowfmoutput/flatbed_quadrilateral_his.nc','waterlevel');

% Define maps
map1        = ['../c030_flatbed_quadrilateral/dflowfmoutput/flatbed_quadrilateral_map.nc'];
map2        = ['../c031_flatbed_quadrilateral_fine/dflowfmoutput/flatbed_quadrilateral_map.nc'];
map3        = ['../c032_flatbed_quadrilateral_finest/dflowfmoutput/flatbed_quadrilateral_map.nc'];

% Read output
y1          = nc_varget(map1,'FlowElem_ycc');
y2          = nc_varget(map2,'FlowElem_ycc');
y3          = nc_varget(map3,'FlowElem_ycc');
h1          = nc_varget(map1,'s1');
h2          = nc_varget(map2,'s1');
h3          = nc_varget(map3,'s1');
h1          = h1(end,:)';
h2          = h2(end,:)';
h3          = h3(end,:)';

% Compute exact solution
B           =  300000;
numseconds  =  23*60*60 + 56*60 + 4.1;
Omega       =  2*pi/numseconds;
g           =  9.81;
phi         =  45;
u           =  0.1;
d           =  500;
ynul        =  150000;
f           =  2.*Omega.*sind(phi);
ha1         = -f.*(y1-ynul)./g.*u;
ha2         = -f.*(y2-ynul)./g.*u;
ha3         = -f.*(y3-ynul)./g.*u;

mean(h1-ha1)
mean(h2-ha2)
mean(h3-ha3)

% Determine L2-norm for Cartesian grids (waterlevel)
dhsqr       = (h1 - ha1).^2;
dhsqrsum    = sum(dhsqr(:));
diffh(1)    = sqrt(dhsqrsum)/sqrt(length(y1));
dhsqr       = (h2 - ha2).^2;
dhsqrsum    = sum(dhsqr(:));
diffh(2)    = sqrt(dhsqrsum)/sqrt(length(y2));
dhsqr       = (h3 - ha3).^2;
dhsqrsum    = sum(dhsqr(:));
diffh(3)    = sqrt(dhsqrsum)/sqrt(length(y3));

% Cartesian grids
dx(1)       = 20e3;
dx(2)       = 10e3;
dx(3)       =  5e3;

% Minus 1 slope
xvgl1       = [1e-5 1e-3];
yvgl1       = xvgl1.^(-1)/1e7;

% Figure settings
ms1         = 25  ;
ms2         =  8  ;
fs          = 20  ;
cl          =  0.5;
lw          =  1  ;

% Make plot (waterlevel)
figure(1);
loglog(1./dx,diffh,'.-','markersize',ms1,'linewidth',lw);    hold on;
loglog(xvgl1,yvgl1,'--','color',[cl cl cl],'linewidth',lw);  hold off;
grid on;
xlim([1e-5 1e-3]);
ylim([1e-4 1e-2]);
text(3.0e-3,4.0e-8,'1st order','color',[cl cl cl],'fontsize',fs);
text(1.0e-2,1.2e-7,'2nd order','color',[cl cl cl],'fontsize',fs);
xlabel('inverse flow link length [m^{-1}]','fontsize',fs);
ylabel('L_2-norm [m]','fontsize',fs);
title('Convergence behavior (water level)','fontsize',fs);
set(gca,'fontsize',fs);

% Print?
if printen == 1;
    print(figure(1),'-dpng','-r300','doc/coriolisstraightconvergence.png');
    close all;
end