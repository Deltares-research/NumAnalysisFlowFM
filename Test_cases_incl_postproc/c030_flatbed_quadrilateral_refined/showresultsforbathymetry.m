function showresultsforbathymetry(shiften);

fclose all;
close all;
clc;

% Configuration
if nargin == 0;
    shiften =  0;
end

% Print?
printen     =  1;

% Define output
if shiften == 0;
    map1    = ['../c030_flatbed_quadrilateral/dflowfmoutput/flatbed_quadrilateral_map.nc'];
    map2    = ['../c034_linearbed_quadrilateral/dflowfmoutput/linearbed_quadrilateral_map.nc'];
    map3    = ['../c035_linearbed_quadrilateral_piecewiselinear/dflowfmoutput/linearbed_quadrilateral_map.nc'];
    map4    = ['../c036_cosine_squares/dflowfmoutput/cosine_squares_map.nc'];
else
    map1    = ['../c030_flatbed_quadrilateral/dflowfmoutput/flatbed_quadrilateral_map.nc'];
    map2    = ['../c040_riemannbndshift_linear/dflowfmoutput/linearbed_quadrilateral_map.nc'];
    map3    = ['../c041_riemannbndshift_piecewiselinear/dflowfmoutput/linearbed_quadrilateral_map.nc'];
    map4    = ['../c042_riemannbndshift_cosine/dflowfmoutput/cosine_squares_map.nc'];
end

% Read output
y1          = nc_varget(map1,'FlowElem_ycc');
y2          = nc_varget(map2,'FlowElem_ycc');
y3          = nc_varget(map3,'FlowElem_ycc');
y4          = nc_varget(map4,'FlowElem_ycc');
h1          = nc_varget(map1,'s1');
h2          = nc_varget(map2,'s1');
h3          = nc_varget(map3,'s1');
h4          = nc_varget(map4,'s1');
h1          = h1(end,:);
h2          = h2(end,:);
h3          = h3(end,:);
h4          = h4(end,:);

yy1 = unique(y1);
yy2 = unique(y2);
yy3 = unique(y3);
yy4 = unique(y4);

% Apply averaging
for i=1:15;
    j1      =  find(y1==yy1(i));
    ya1(i)  =  yy1(i);
    ha1(i)  =  mean(h1(j1));
    j2      =  find(y2==yy2(i));
    ya2(i)  =  yy2(i);
    ha2(i)  =  mean(h2(j2));
    j3      =  find(y3==yy3(i));
    ya3(i)  =  yy3(i);
    ha3(i)  =  mean(h3(j3));
    j4      =  find(y4==yy4(i));
    ya4(i)  =  yy4(i);
    ha4(i)  =  mean(h4(j4));
end

% Exact solution
B           =  300000;
numseconds  =  23*60*60 + 56*60 + 4.1;
Omega       =  2*pi/numseconds;
g           =  9.81;
phi         =  45;
u           =  0.1;
d           =  500;
ynul        =  150000;
f           =  2.*Omega.*sind(phi);
y           = [0:1000:300000];
zetafinal   = -f.*(y-ynul)./g.*u;

% Figure settings
fs          = 18;
lw          =  2;
ms          =  20;
cl          =  0.7;

% Plot outcomes
plot(y /1000,zetafinal,'-' ,'linewidth' ,lw,'color',[cl cl cl]); hold on;
plot(ya1/1000,       ha1,'.k','markersize',ms);
plot(ya2/1000,       ha2,'.r','markersize',ms);
plot(ya3/1000,       ha3,'.b','markersize',ms);
plot(ya4/1000,       ha4,'.m','markersize',ms); hold off;
grid on;

% Legend
legend('exact','flat','linear','piecewise','cosine','location','northeast');

% Beautify
xlim([0 B/1000]);
ylim([-0.2 0.2]);
xlabel('lateral orientation [km]','fontsize',fs);
ylabel('vertical orientation [m w.r.t. ref]','fontsize',fs);
set(gca,'fontsize',fs);

% Printen
if printen == 1 & shiften == 0;
    print(gcf,'-dpng','-r300','doc/resultsbathymetry.png');
    close all;
end
if printen == 1 & shiften == 1;
    print(gcf,'-dpng','-r300','doc/resultsbathymetryshift.png');
    close all;
end