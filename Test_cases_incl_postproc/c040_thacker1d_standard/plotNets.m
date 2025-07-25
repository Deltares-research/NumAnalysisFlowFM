function plotNets;

fclose all; clear all; close all; clc;

% Geef bestandsnaam eenmalig
bestand        = ['dflowfmoutput\planar1d_map.nc'];

% Lees netnodes
xn             = nc_varget(bestand,'NetNode_x');
yn             = nc_varget(bestand,'NetNode_y');
netlinks       = nc_varget(bestand,'NetLink');

% Bepaal netlinks
xnn            = xn(netlinks);
ynn            = yn(netlinks);

% Geef figuur settings
ms             = 4;
fs             = 12;
lw             = 0.5;
c1             = 0.3;
c2             = 0.8;
c3             = 0.1;

% Maak plotje van rooster
figure(1)
plot(xn,yn,'.','markersize',ms,'color','k'); hold on;
daspect([1 1 1]);
line(xnn',ynn','color','k','linewidth',lw);
hold off
axis off

% Set bounds
print('-dpng',['doc/grid.png']);
close all;