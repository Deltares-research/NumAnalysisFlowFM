function showgridtypeeffect;

fclose all;
close all;
clc;

% Define maps
map1        = ['../c031_flatbed_quadrilateral_fine/dflowfmoutput/flatbed_quadrilateral_map.nc'];
map2        = ['../c033_flatbed_triangular/dflowfmoutput/flatbed_triangular_map.nc'];
map3        = ['../c036_cosine_squares/dflowfmoutput/cosine_squares_map.nc'];
map4        = ['../c037_cosine_triangles/dflowfmoutput/cosine_triangles_map.nc'];

% Read output
y1          = nc_varget(map1,'FlowElem_ycc');
y2          = nc_varget(map2,'FlowElem_ycc');
y3          = nc_varget(map3,'FlowElem_ycc');
y4          = nc_varget(map4,'FlowElem_ycc');
h1          = nc_varget(map1,'s1');
h2          = nc_varget(map2,'s1');
h3          = nc_varget(map3,'s1');
h4          = nc_varget(map4,'s1');
h1          = h1(end,:)';
h2          = h2(end,:)';
h3          = h3(end,:)';
h4          = h4(end,:)';

% Compute exact solution
B           =  300000;
numseconds  =  23*60*60 + 56*60 + 4.1;
Omega       =  2*pi/numseconds;
g           =  9.81;
phi         =  45;
u           =  0.1;
d           =  500;
ynul        =  150000;
f           =  2.*Omega.*sind(phi);
ha1         = -f.*(y1-ynul)./g.*u;
ha2         = -f.*(y2-ynul)./g.*u;
ha3         = -f.*(y3-ynul)./g.*u;
ha4         = -f.*(y4-ynul)./g.*u;

% Determine L2-norm for Cartesian grids (waterlevel)
dhsqr       = (h1 - ha1).^2;
dhsqrsum    = sum(dhsqr(:));
diffh(1)    = sqrt(dhsqrsum)/sqrt(length(y1));
dhsqr       = (h2 - ha2).^2;
dhsqrsum    = sum(dhsqr(:));
diffh(2)    = sqrt(dhsqrsum)/sqrt(length(y2));
dhsqr       = (h3 - ha3).^2;
dhsqrsum    = sum(dhsqr(:));
diffh(3)    = sqrt(dhsqrsum)/sqrt(length(y3));
dhsqr       = (h4 - ha4).^2;
dhsqrsum    = sum(dhsqr(:));
diffh(4)    = sqrt(dhsqrsum)/sqrt(length(y4));

diffh