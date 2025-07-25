function [x,h] = belanger(Q,C,L,B,ib,hb,dxfm);

% Number of increments
N         = L/dxfm*1000;

% Set constant
g         = 9.81;

% Compute he and hg
hg        = (Q.^2./B.^2./g).^(1/3);
he        = (Q.^2./B.^2./C.^2./ib).^(1/3);

% Define grid
dx        = L/N;
x         = [0:dx:L+dxfm];
sx        = length(x);

% Compute surface level
h(sx)     = hb;
for i = (sx-1):-1:1;
    rhs   = ib.*(h(i+1).^3-he.^3)./(h(i+1).^3-hg.^3);
    h(i)  = h(i+1) - dx.*rhs;
end

% Figure settings
fs        = 14;
lw        = 1;