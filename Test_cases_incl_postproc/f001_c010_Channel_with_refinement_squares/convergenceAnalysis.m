function convergenceAnalysis;

% Clean up
clc;
close all;
warning off all;
format long;

% Input parameters
L         = 10000                  ;
Q         =  2500                  ;
B         =   500                  ;
C         =    65                  ;
deltab    =     1                  ;
ib        = deltab./L              ;
deq       = (Q/B/C/sqrt(ib)).^(2/3);
ueq       = (Q/B/deq)              ;

% Read water levels (=> 'hu'  at flow links  )
[tabel]   = computeBedlevelAtFlowlink('dflowfmoutput\simplechannel_map.nc', 4);

% Read velocities   (=> 'ucx' at cell centers)
yf        = nc_varget('dflowfmoutput/simplechannel_map.nc','FlowElem_ycc');
usim      = nc_varget('dflowfmoutput/simplechannel_map.nc','ucx');
vsim      = nc_varget('dflowfmoutput/simplechannel_map.nc','ucy');
usim      = usim(end,:)';
vsim      = vsim(end,:)';

% Find grids (for flow links)
grid1     = find(tabel(:,11)>=    0 & tabel(:,11)<=  500); 
grid2     = find(tabel(:,11)>= 1000 & tabel(:,11)<= 1500); 
grid3     = find(tabel(:,11)>= 2000 & tabel(:,11)<= 2500); 
grid4     = find(tabel(:,11)>= 3000 & tabel(:,11)<= 3500); 
grid5     = find(tabel(:,11)>= 4000 & tabel(:,11)<= 4500); 

% Determine L2-norm for Cartesian grids (waterlevel)
dhsqr     = (tabel(grid1,12) - deq).^2;
dhsqrsum  = sum(dhsqr(:));
diffc(1)  = sqrt(dhsqrsum)/sqrt(length(grid1));
dhsqr     = (tabel(grid2,12) - deq).^2;
dhsqrsum  = sum(dhsqr(:));
diffc(2)  = sqrt(dhsqrsum)/sqrt(length(grid2));
dhsqr     = (tabel(grid3,12) - deq).^2;
dhsqrsum  = sum(dhsqr(:));
diffc(3)  = sqrt(dhsqrsum)/sqrt(length(grid3));
dhsqr     = (tabel(grid4,12) - deq).^2;
dhsqrsum  = sum(dhsqr(:));
diffc(4)  = sqrt(dhsqrsum)/sqrt(length(grid4));
dhsqr     = (tabel(grid5,12) - deq).^2;
dhsqrsum  = sum(dhsqr(:));
diffc(5)  = sqrt(dhsqrsum)/sqrt(length(grid5));

% Find grids (for cell centers)
grid1     = find(yf>=    0 & yf<=  500); 
grid2     = find(yf>= 1000 & yf<= 1500); 
grid3     = find(yf>= 2000 & yf<= 2500); 
grid4     = find(yf>= 3000 & yf<= 3500); 
grid5     = find(yf>= 4000 & yf<= 4500);

% Determine L2-norm for Cartesian grids (velocity)
dusqr     = (usim(grid1) - ueq).^2 + (vsim(grid1) - 0).^2;
dusqrsum  = sum(dusqr(:));
diffcu(1) = sqrt(dusqrsum)/sqrt(length(grid1));
dusqr     = (usim(grid2) - ueq).^2 + (vsim(grid2) - 0).^2;
dusqrsum  = sum(dusqr(:));
diffcu(2) = sqrt(dusqrsum)/sqrt(length(grid2));
dusqr     = (usim(grid3) - ueq).^2 + (vsim(grid3) - 0).^2;
dusqrsum  = sum(dusqr(:));
diffcu(3) = sqrt(dusqrsum)/sqrt(length(grid3));
dusqr     = (usim(grid4) - ueq).^2 + (vsim(grid4) - 0).^2;
dusqrsum  = sum(dusqr(:));
diffcu(4) = sqrt(dusqrsum)/sqrt(length(grid4));
dusqr     = (usim(grid5) - ueq).^2 + (vsim(grid5) - 0).^2;
dusqrsum  = sum(dusqr(:));
diffcu(5) = sqrt(dusqrsum)/sqrt(length(grid5));

% Cartesian grids
dc(1)     = 500.00;
dc(2)     = 250.00;
dc(3)     = 125.00;
dc(4)     =  62.50;
dc(5)     =  31.25;

% Figure settings
ms1       = 25  ;
ms2       =  8  ;
fs        = 20  ;
cl        =  0.5;
lw        =  1  ;

% Minus 1 slope
xvgl1     = [1e-4 1e0];
yvgl1     = xvgl1.^(-1)/1e10;

% Minus 2 slope
xvgl2     = [1e-4 1e0];
yvgl2     = xvgl2.^(-2)/1e11;

% Make plot (waterlevel)
figure(1);
loglog(1./dc,diffc,'.-','markersize',ms1,'linewidth',lw);    hold on;
loglog(xvgl1,yvgl1,'--','color',[cl cl cl],'linewidth',lw);  
loglog(xvgl2,yvgl2,'--','color',[cl cl cl],'linewidth',lw);  hold off;
grid on;
xlim([1e-3 1e-1]);
ylim([1e-9 1e-6]);
text(3.0e-3,4.0e-8,'1st order','color',[cl cl cl],'fontsize',fs);
text(1.0e-2,1.2e-7,'2nd order','color',[cl cl cl],'fontsize',fs);
legend('Cartesian grids','location','southeast');
xlabel('inverse flow link length [m^{-1}]','fontsize',fs);
ylabel('L_2-norm [m]','fontsize',fs);
title('Convergence behavior (water depth)','fontsize',fs);
set(gca,'fontsize',fs);

% Minus 1 slope
xvgl1     = [1e-4 1e0];
yvgl1     = xvgl1.^(-1)/1e10;

% Minus 2 slope
xvgl2     = [1e-4 1e0];
yvgl2     = xvgl2.^(-2)/1e11;

% Make plot (velocity)
figure(2);
loglog(1./dc,diffcu,'.-','markersize',ms1,'linewidth',lw);   hold on;
loglog(xvgl1,yvgl1,'--','color',[cl cl cl],'linewidth',lw);  
loglog(xvgl2,yvgl2,'--','color',[cl cl cl],'linewidth',lw);  hold off;
grid on;
xlim([1e-3 1e-1]);
ylim([1e-9 1e-6]);
text(3.0e-3,4.0e-8,'1st order','color',[cl cl cl],'fontsize',fs);
text(1.0e-2,1.2e-7,'2nd order','color',[cl cl cl],'fontsize',fs);
legend('Cartesian grids','location','northeast');
xlabel('inverse flow link length [m^{-1}]','fontsize',fs);
ylabel('L_2-norm [m]','fontsize',fs);
title('Convergence behavior (flow velocity)','fontsize',fs);
set(gca,'fontsize',fs);

% Print
print(figure(1),'-dpng','-r300','doc/hchannelconvergence.png');
print(figure(2),'-dpng','-r300','doc/uchannelconvergence.png');
close all;