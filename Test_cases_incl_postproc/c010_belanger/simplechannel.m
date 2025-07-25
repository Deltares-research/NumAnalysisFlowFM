clc;
close all;
warning off all;
format long;

% Input
C         = 65;
Q         = 600;
L         = 100000;
B         = 20;
ib        = 1e-4;
dx        = 500;

% Boundary condition
hb        =  0;
db        = hb + ib.*(L + dx);

% Analytical solution
[xb,heq]  = belanger(Q,C,L,B,ib,db,dx);

% Read waterlevel
xf        = nc_varget('dflowfmoutput/simplechannel_map.nc','FlowElem_xcc');
hsim      = nc_varget('dflowfmoutput/simplechannel_map.nc','s1'          );
hsim      = hsim(end,:)';

% Locations for water levels actually at netnodes
xn        = xf + dx/2;

% Compute waterdepth
ds        = -xn.*ib;
hsim      = hsim -ds;

% Find computation points in analytical solution
for i = 1:length(xn);
    j          = find(xb==xn(i));
    heqs(i,1)  = heq(j);
end
ddiff     = sqrt(mean((heqs - hsim).^2));

% Figure settings
ms1       = 15  ;
ms2       =  8  ;
ms3       = 25  ;
fs        = 16  ;
cl        =  0.5;
lw        =  1  ;

% Make plot
figure(1);
plot(xn/1e3 ,heqs ,'k-' ,'linewidth',lw,'markersize',ms2); hold on;
plot(xn/1e3 ,hsim ,'r-' ,'linewidth',lw,'markersize',ms2); 
plot((L + dx)/1e3 ,db   ,'k.' ,'linewidth',lw,'markersize',ms3); hold off;
grid on;

% Set labels and legend
xlim([0 102]);
legend('semi-analytical solution','dflowfm','boundary condition','location','southwest');
xlabel('distance from inflow boundary [km]','fontsize',fs);
ylabel('water depth [m]','fontsize',fs);
text(16.5,11.7,['rms difference = ',num2str(ddiff,'%1.2e'),' m'],'color',[cl cl cl],'fontsize',fs);
title('Waterdepth along the channel','fontsize',fs);
set(gca,'fontsize',fs);

% Print
print(figure(1),'-dpng','-r300','doc/waterdepth.png');
close all;