function showbathymetry;

fclose all;
close all;
clc;

% Print?
printen        = 1;

% Input
d1             =     50;
d2             =    500;
B              = 300000;

% Set grid in y-direction
y              = [0:1000:300000];

% Define four bathymetries
z1             = -d2.*ones(1,length(y));
z2             = -(d2-d1)./B.*y - 50;
z3             = -abs((d1-d2)/(B-2.*40000).*(y-40000)-d1);
z4             = -abs((d2-d1)./2.*cos(2.*pi.*(y)./2./B) - (d2-d1)./2 - d1);

% Modify z3
z3(y<= 40000)  = -d1;
z3(y>=260000)  = -d2;

% Solid parameters
numseconds     =  23*60*60 + 56*60 + 4.1;
Omega          =  2*pi/numseconds;
g              =  9.81;

% Input
phi            =  45;
u              =  0.1;
d              =  500;
ynul           =  150000;

% Preprocessing
f              =  2.*Omega.*sind(phi);

% Surface elevation
zetafinal      = -f.*(y-ynul)./g.*u;

% Figure settings
fs             = 18;
lw             =  1.5;

% Make plot
plot(y/1000,z1,'k','linewidth',lw);  hold on;
plot(y/1000,z2,'r','linewidth',lw); 
plot(y/1000,z3,'b','linewidth',lw);
plot(y/1000,z4,'m','linewidth',lw); 
plot(y/1000,zetafinal,'c'); hold off;
grid on;

% Legend
legend('flat bed','linearly varying','piecewise linearly varying','cosine','location','northeast');

% Beautify
xlim([0 B/1000]);
ylim([-600 100]);
xlabel('lateral orientation [km]','fontsize',fs);
ylabel('vertical orientation [m]','fontsize',fs);
set(gca,'fontsize',fs);

% Print?
if printen == 1;
    print(gcf,'-dpng','-r300','doc/coriolisstraightconfig.png');
    close all;
end