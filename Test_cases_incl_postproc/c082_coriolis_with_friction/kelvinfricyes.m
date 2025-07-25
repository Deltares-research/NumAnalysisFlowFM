fclose all;
clear all;
close all;
clc;

% Geometry
B            =  550;
L            = 1800;

% Load analytical solution
load('kelvinfricyes.mat');
xana         = xana/1e3;
yana         = yana/1e3;
[xana,yana]  = meshgrid(xana,yana);
yana         = yana + B/2;

% Read the map-file
map          = ['dflowfmoutput/kelvinfricyes_map.nc'];
h            = nc_varget(map,'s1');
for i=1:size(h,2);
    hmax(i)  = max(h(:,i));
end

% Load the geometry (flow nodes)
x            = nc_varget(map,'FlowElem_xcc');
y            = nc_varget(map,'FlowElem_ycc');

% Plot analytical solution
pcolor(xana,yana,hana);
axis off;
shading flat;
daspect([1 1 1]);

% Read solution
G            = dflowfm.readNet(map);
D            = hmax;

% Plot first solution
hold on;
fac          = 1.3;
xsol         = G.cor.x/1e3;
ysol         = G.cor.y/1e3 + fac.*B;
h            = trisurfcorcen(G.tri,xsol,ysol,D(G.map3),D(G.map3));
set(h,'edgeColor','none');
hold off;

% Beautify
xlim([0 L]);
ylim([0 max(ysol)]);
caxis([0 1.6])
daspect([1 1 1]);

% Print
print(figure(1),'-dpng','-r300','doc/kelvinfricyes.png');
close all;