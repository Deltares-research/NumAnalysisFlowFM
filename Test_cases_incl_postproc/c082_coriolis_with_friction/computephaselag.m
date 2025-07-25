clear all; 
close all;
fclose all;
clc;

g         =   9.81;

% Geometry
L         = 1800e3;
B         =  550e3;
H         =     80;

% Parameters k en omega
T         = 745*60;
omega     = 2*pi/T;
k         = omega./sqrt(g*H);

% Input for formules
x         = L;
t         = [0:T/3e3:T];

sigA      = cos(k.*x).*cos(omega.*t) - sin(k.*x).*sin(omega.*t);
sigB      = cos(k.*x).*cos(omega.*t) + sin(k.*x).*sin(omega.*t);

plot(t/T,sigA,'b'); hold on;
plot(t/T,sigB,'r'); hold off;
grid on;
xlim([0 1]);

maxdifft  = t(sigA==max(sigA(:))) - t(sigB==max(sigB(:)));
mindifft  = t(sigA==min(sigA(:))) - t(sigB==min(sigB(:)));

difft1    = maxdifft/T;
difft2    = mindifft/T;

fasediff1 = difft1*360
fasediff2 = difft2*360