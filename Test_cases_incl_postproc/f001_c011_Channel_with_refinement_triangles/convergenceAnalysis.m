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
xf        = nc_varget('dflowfmoutput/simplechannel_map.nc','FlowElem_xcc');
yf        = nc_varget('dflowfmoutput/simplechannel_map.nc','FlowElem_ycc');
usim      = nc_varget('dflowfmoutput/simplechannel_map.nc','ucx');
vsim      = nc_varget('dflowfmoutput/simplechannel_map.nc','ucy');
usim      = usim(end,:)';
vsim      = vsim(end,:)';

% Find grids (for flow links)
grid6     = find(tabel(:,11)>= 5000 & tabel(:,11)<= 5500);
grid7     = find(tabel(:,11)>= 6000 & tabel(:,11)<= 6500);
grid8     = find(tabel(:,11)>= 7000 & tabel(:,11)<= 7500);

% Consider medium part
grid6med  = find(tabel(:,10)>=  3000 & tabel(:,10)<=  7000 & tabel(:,11)>= 5000 & tabel(:,11)<= 5500); 
grid7med  = find(tabel(:,10)>=  3000 & tabel(:,10)<=  7000 & tabel(:,11)>= 6000 & tabel(:,11)<= 6500); 
grid8med  = find(tabel(:,10)>=  3000 & tabel(:,10)<=  7000 & tabel(:,11)>= 7000 & tabel(:,11)<= 7500); 
g6        = mean(tabel(grid6med,12));
g7        = mean(tabel(grid7med,12));
g8        = mean(tabel(grid8med,12));
display(['Deviations for the medium part of the domain:']);
display(['   coarse grid: ',num2str(g6,'%3.8f'),' - ',num2str(deq,'%3.8f'),' =  ',num2str(g6-deq,'%3.8f')]);
display(['   medium grid: ',num2str(g7,'%3.8f'),' - ',num2str(deq,'%3.8f'),' = ',num2str(g7-deq,'%3.8f')]);
display(['   fine   grid: ',num2str(g8,'%3.8f'),' - ',num2str(deq,'%3.8f'),' =  ',num2str(g8-deq,'%3.8f')]);

% Determine L2-norm for triangular grids (waterlevel): with boundaries
dhsqr     = (tabel(grid6,12) - deq).^2;
dhsqrsum  = sum(dhsqr(:));
difft(1)  = sqrt(dhsqrsum)/length(grid6);
dhsqr     = (tabel(grid7,12) - deq).^2;
dhsqrsum  = sum(dhsqr(:));
difft(2)  = sqrt(dhsqrsum)/length(grid7);
dhsqr     = (tabel(grid8,12) - deq).^2;
dhsqrsum  = sum(dhsqr(:));
difft(3)  = sqrt(dhsqrsum)/length(grid8);

% Find grids (for cell centers)
grid6     = find(yf>= 5000 & yf<= 5500); 
grid7     = find(yf>= 6000 & yf<= 6500); 
grid8     = find(yf>= 7000 & yf<= 7500);

% Consider medium part
grid6med  = find(xf>=  3000 & xf<=  7000 & yf>= 5000 & yf<= 5500); 
grid7med  = find(xf>=  3000 & xf<=  7000 & yf>= 6000 & yf<= 6500); 
grid8med  = find(xf>=  3000 & xf<=  7000 & yf>= 7000 & yf<= 7500); 

% Determine L2-norm for triangular grids (velocity): with boundaries
dusqr     = (usim(grid6) - ueq).^2 + (vsim(grid6) - 0).^2;
dusqrsum  = sum(dusqr(:));
difftu(1) = sqrt(dusqrsum)/length(grid6);
dusqr     = (usim(grid7) - ueq).^2 + (vsim(grid7) - 0).^2;
dusqrsum  = sum(dusqr(:));
difftu(2) = sqrt(dusqrsum)/length(grid7);
dusqr     = (usim(grid8) - ueq).^2 + (vsim(grid8) - 0).^2;
dusqrsum  = sum(dusqr(:));
difftu(3) = sqrt(dusqrsum)/length(grid8);

% Triangular grids
dt(1)     = 500/3/1;
dt(2)     = 500/3/2;
dt(3)     = 500/3/4;

% Figure settings
ms1       = 25  ;
ms2       =  8  ;
fs        = 20  ;
cl        =  0.5;
lw        =  1  ;

% Minus 1 slope
xvgl1     = [1e-4 1e0];
yvgl1     = xvgl1.^(-1)/1e5;

% Minus 2 slope
xvgl2     = [1e-4 1e0];
yvgl2     = xvgl2.^(-2)/1e6;

% Make plot (waterlevel)
figure(1);
loglog(1./dt,difft,'o-','markersize',ms2,'linewidth',lw);    hold on;
loglog(xvgl1,yvgl1,'--','color',[cl cl cl],'linewidth',lw);  
loglog(xvgl2,yvgl2,'--','color',[cl cl cl],'linewidth',lw);  hold off;
grid on;
xlim([1e-3 1e-1]);
ylim([1e-6 1e-1]);
text(1.4e-3,9.0e-3,'1st order','color',[cl cl cl],'fontsize',fs);
text(1.0e-2,1.2e-2,'2nd order','color',[cl cl cl],'fontsize',fs);
legend('Triangular grids','location','southwest');
xlabel('inverse flow link length [m^{-1}]','fontsize',fs);
ylabel('L_2-norm [m]','fontsize',fs);
title('Convergence behavior (water depth)','fontsize',fs);
set(gca,'fontsize',fs);

% Make plot (velocity)
figure(2);  
loglog(1./dt,difftu,'o-','markersize',ms2,'linewidth',lw);   hold on;
loglog(xvgl1,yvgl1,'--','color',[cl cl cl],'linewidth',lw);  
loglog(xvgl2,yvgl2,'--','color',[cl cl cl],'linewidth',lw);  hold off;
grid on;
xlim([1e-3 1e-1]);
ylim([1e-6 1e-1]);
text(1.4e-3,9.0e-3,'1st order','color',[cl cl cl],'fontsize',fs);
text(1.0e-2,1.2e-2,'2nd order','color',[cl cl cl],'fontsize',fs);
legend('Triangular grids','location','southwest');
xlabel('inverse flow link length [m^{-1}]','fontsize',fs);
ylabel('L_2-norm [m]','fontsize',fs);
title('Convergence behavior (velocity)','fontsize',fs);
set(gca,'fontsize',fs);

% Print
print(figure(1),'-dpng','-r300','doc/hchannelconvergence.png');
print(figure(2),'-dpng','-r300','doc/uchannelconvergence.png');
close all;