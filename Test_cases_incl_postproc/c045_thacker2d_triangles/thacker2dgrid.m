function sobeygrid;

fclose all;
close all;
clear all;
clc;


% Print?
if nargin == 0;
    printen = 1;
end

% Define output file
bestand        = ['radial2dtri_net.nc'];

% Read grid
xn             = nc_varget(bestand,'NetNode_x');
yn             = nc_varget(bestand,'NetNode_y');
netlinks       = nc_varget(bestand,'NetLink');
xnn            = xn(netlinks);
ynn            = yn(netlinks);

% Figure settings
ms1       = 4;
ms2       = 20;
fs        = 18;
lw        = 0.5;
cl        = 0.5;

% Make plot of grid
figure(1)
plot(xn,yn,'.','markersize',ms1,'color',[cl cl cl]); hold on;
line(xnn',ynn','color',[cl cl cl],'linewidth',lw);
daspect([1 1 1]);

% Confine domain
R         = 12e4;
xmin      =  6e4;
ymin      =  6e4;
xmax      = 11e4;
ymax      = 10e4;
xlim([xmin xmax]);
ylim([ymin ymax]);

% Beautify
set(gca,'fontsize',fs);
xc        = [xmin:100:xmax];
yc        = sqrt(R.^2 - xc.^2);
plot(xc,yc,'k','linewidth',lw*3);
xlabel('x-coordinate [m]','fontsize',fs);
ylabel('y-coordinate [m]','fontsize',fs);

% Print
if printen == 1;
    print(figure(1),'-dpng','-r300','doc/thacker2dgridtriangles.png');
    close all;
end