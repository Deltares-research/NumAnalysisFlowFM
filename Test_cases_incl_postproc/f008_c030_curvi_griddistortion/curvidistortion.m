fclose all;
clear all;
close all;
clc;

% Read the map-file
map          = ['dflowfmoutput/curvidistortion_map.nc'];
hread        = nc_varget(map,'s1');
hsol         = hread(end,:)';

% Load the geometry (flow nodes)
x            = nc_varget(map,'FlowElem_xcc');
y            = nc_varget(map,'FlowElem_ycc');

% Read solution
G            = dflowfm.readNet(map);

% Read depth file
dep          = load('depth.xyz');

% Plot first solution
hold on;
xsol         = G.cor.x;
ysol         = G.cor.y;
dsol         = dep(1,3) + (dep(3,3)-dep(1,3))./(dep(3,1)-dep(1,1)).*x;
D            = hsol - dsol;
hfig         = trisurfcorcen(G.tri,xsol,ysol,D(G.map3));
set(hfig,'edgeColor','none');
hold off;

% Figure settings
fs           = 14;

% Beautify
xlim([min(xsol(:)) max(xsol(:))]);
ylim([0 540]);
box on;
daspect([1 1 1]);
caxis([min(D(:)) max(D(:))]);
colorbar('ytick',[min(D(:)) 2.01 2.02 2.03 2.037 max(D(:))]);
set(gca,'ytick',[0 240 300 540]);
set(gca,'yticklabel',[{'0'} {'240'} {'0'} {'240'}]);
xlabel('longitudinal distance [m]','fontsize',fs);
ylabel('lateral distance [m]','fontsize',fs);
title('water depth [m]','fontsize',fs);
set(gca,'fontsize',fs);

% Print
print(figure(1),'-dpng','-r300','doc/distortedgridresult.png');
close all;