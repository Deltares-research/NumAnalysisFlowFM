function showgrids;

fclose all;
close all;
clc;

% Printen?
printen     = 1;

% Gridnames
grid1       = ['../c036_cosine_squares/cosine_squares_net.nc'];
grid2       = ['../c037_cosine_triangles/cosine_triangles_net.nc'];

% Read grid1
xn1         = nc_varget(grid1,'NetNode_x');
yn1         = nc_varget(grid1,'NetNode_y');
zn1         = nc_varget(grid1,'NetNode_z');
netlinks1   = nc_varget(grid1,'NetLink');
xnn1        = xn1(netlinks1);
ynn1        = yn1(netlinks1);

% Read grid2
xn2         = nc_varget(grid2,'NetNode_x');
yn2         = nc_varget(grid2,'NetNode_y');
zn2         = nc_varget(grid2,'NetNode_z');
netlinks2   = nc_varget(grid2,'NetLink');
xnn2        = xn2(netlinks2);
ynn2        = yn2(netlinks2);

% Figure settings
ms1         = 10;
ms2         = 50;
fs          = 18;
lw          = 1;
cl          = 0.7;

% Make plot of grid1
figure(1)
line(xnn1',ynn1','color',[cl cl cl],'linewidth',lw); hold on;
plot(xn1,yn1,'.','markersize',ms1,'color',[cl cl cl]); 
scatter3(xn1,yn1,zn1,ms2,zn1,'filled');
colormap jet
caxis([-500 -50]);
h = colorbar;
set(h,'ytick',[-500 -50]);
daspect([1 1 1]);
xlim([0e5 5e5]);
ylim([0e5 3e5]);
hold off
xlabel('longitudinal orientation [m]','fontsize',fs);
ylabel('lateral orientation [m]','fontsize',fs);
set(gca,'fontsize',fs);

% Make plot of grid2
figure(2)
line(xnn2',ynn2','color',[cl cl cl],'linewidth',lw); hold on;
plot(xn2,yn2,'.','markersize',ms1,'color',[cl cl cl]); 
scatter3(xn2,yn2,zn2,ms2,zn2,'filled');
colormap jet
caxis([-500 -50]);
h = colorbar;
set(h,'ytick',[-500 -50]);
daspect([1 1 1]);
xlim([0e5 5e5]);
ylim([0e5 3e5]);
hold off
xlabel('longitudinal orientation [m]','fontsize',fs);
ylabel('lateral orientation [m]','fontsize',fs);
set(gca,'fontsize',fs);

% Print?
if printen == 1;
    print(figure(1),'-dpng','-r300','doc/coriolisstraightsquares.png');
    print(figure(2),'-dpng','-r300','doc/coriolisstraighttriangles.png');
    close all;
end