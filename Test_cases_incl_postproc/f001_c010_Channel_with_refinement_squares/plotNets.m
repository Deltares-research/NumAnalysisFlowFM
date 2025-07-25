function plotNets;

fclose all; clear all; close all; clc;

%% Geef bestandsnaam eenmalig
bestand        = ['dflowfmoutput\simplechannel_map.nc'];

%% Lees netnodes
xn             = nc_varget(bestand,'NetNode_x');
yn             = nc_varget(bestand,'NetNode_y');
netlinks       = nc_varget(bestand,'NetLink');

%% Lees flownodes
xf             = nc_varget(bestand,'FlowElem_xcc');
yf             = nc_varget(bestand,'FlowElem_ycc');
flowlinks      = nc_varget(bestand,'FlowLink');
nflowlinks     = length(xf);
flowlinks(find(flowlinks>nflowlinks),:) = [];

%% Bepaal netlinks
xnn            = xn(netlinks);
ynn            = yn(netlinks);

%% Bepaal flowlinks
xff            = xf(flowlinks);
yff            = yf(flowlinks);

%% Geef figuur settings
ms             = 4;
fs             = 12;
lw             = 1;
c1             = 0.3;
c2             = 0.8;
c3             = 0.1;

%% Maak plotje van rooster
figure(1)
plot(xn,yn,'.','markersize',ms,'color','k'); hold on;
% plot(xf,yf,'.','markersize',ms,'color','b'); hold off;
daspect([1 1 1]);
% hold on;
line(xnn',ynn','color','k','linewidth',lw);
% line(xff',yff','color','b','linewidth',lw);
hold off
axis off

% %% Set bounds
% print(figure(1),'-dpng','-r300','doc/multiplesimplechannels.png');
% close all;
