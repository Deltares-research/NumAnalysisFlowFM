function thacker1d_planar;

% Basic stuff
close all; 
clc;

% Choice
choice = 0;

% Simulation
if choice == 1;
    sim    = ['thacker1d'];
else
    sim    = ['planar1d'];
end

% Read stuff
xf     = nc_varget(['dflowfmoutput/',sim,'_map.nc'],'FlowElem_xcc');
yf     = nc_varget(['dflowfmoutput/',sim,'_map.nc'],'FlowElem_ycc');
s1     = nc_varget(['dflowfmoutput/',sim,'_map.nc'],'s1');
time   = nc_varget(['dflowfmoutput/',sim,'_map.nc'],'time');

% Grid 1
j      = find(yf<10);
x1     = xf(j);
y1     = yf(j);
h1     = s1(:,j);

% Grid 2
j      = find(yf>10 & yf<20);
x2     = xf(j);
y2     = yf(j);
h2     = s1(:,j);

% Grid 3
j      = find(yf>20);
x3     = xf(j);
y3     = yf(j);
h3     = s1(:,j);

% Basic parameters
g      = 9.81;
lw     = 2;
ms     = 6;
fs     = 18;

% Parameters
L      =  300;
r0     =  120.0;
dep    =  10;

% Derived quantities
omega  = sqrt(2.*g.*dep)./r0;
per    = 2.*pi./omega;

% Geometry
dL     = L/10000;
x      = [-L/2:dL:L/2];
xf     = xf - L/2;
x1     = x1 - L/2;
x2     = x2 - L/2;
x3     = x3 - L/2;
z      = -dep.*(1 - x.^2./r0.^2);
z1     = -dep.*(1 - x1.^2./r0.^2);
z2     = -dep.*(1 - x2.^2./r0.^2);
z3     = -dep.*(1 - x3.^2./r0.^2);
zf     = -dep.*(1 - xf.^2./r0.^2);

% Second input
eta    = 0.23;
psi    = eta.*r0;
samp   = eta.*dep./r0;

% Figure settings
figure(1);
set(gca,'fontsize',fs);

% Start loop
tim    = 1;
sper   = 10;
nper   = 5;
dt     = 2;
for t  = 0:dt:260;

    % Waterlevel
    h      = samp.*(2.*x.*cos(omega.*t) - psi.*(cos(omega.*t)).^2);
    hf     = s1(tim,:);
    h1t    = h1(tim,:);
    j      = find(h1t<=z1');
    h1t(j) = NaN;
    h2t    = h2(tim,:);
    j      = find(h2t<=z2');
    h2t(j) = NaN;
    h3t    = h3(tim,:);
    j      = find(h3t<=z3');
    h3t(j) = NaN;
    j      = find(hf<=zf');
    hf(j)  = NaN;
    j      = find(h<z);
    h(j)   = NaN;
    
    % Make plot
    plot(x ,z  ,'k','linewidth',lw); hold on;
    plot(x ,h  ,'b','linewidth',lw);
    plot(x1,h1t,'r-','linewidth',lw);
    plot(x2,h2t,'r--','linewidth',lw,'markersize',ms); 
    plot(x3,h3t,'r.','linewidth',lw); hold off;
    grid on;
    xlim([-L/2 L/2]);
    title(['Water level after ',num2str(t/per,'%2.2f'),' periods']);
    drawnow;
    pause(0.0);
    
    tim    = tim + 1;

end

% Labels
xlabel('Distance to deepest point in basin [m]','fontsize',fs);
ylabel('Vertical orientation [m]','fontsize',fs);
legend('basin','exact solution','2D grid','1D grid','1D network','location','southeast');

if choice == 0;
    
    % Read waterlevel from his file
    time   = nc_varget(['dflowfmoutput/',sim,'_his.nc'],'time');
    s1     = nc_varget(['dflowfmoutput/',sim,'_his.nc'],'waterlevel');
    u1     = nc_varget(['dflowfmoutput/',sim,'_his.nc'],'x_velocity');
    
    % Generate analytical solution
    x      = 211.5 - L/2;         % x from observation point file
    h      = samp.*(2.*x.*cos(omega.*time) - psi.*(cos(omega.*time)).^2);
    u      = -psi*omega*sin(omega*time);
    
    % Plot the water level
    figure(2);
    plot(time,h ,'b','linewidth',lw); hold on;
    plot(time,s1,'r','linewidth',lw); hold off;
    grid on;
    legend('exact solution','1D grid','location','west');
    set(gca,'fontsize',fs);
    xlabel('Time [seconds]','fontsize',fs);
    ylabel('Surface elevation [m]','fontsize',fs);
    
    % Plot the velocity
    figure(3);
    plot(time,u ,'b','linewidth',lw); hold on;
    plot(time,u1,'r','linewidth',lw); hold off;
    grid on;
    legend('exact solution','1D grid','location','west');
    set(gca,'fontsize',fs);
    xlabel('Time [seconds]','fontsize',fs);
    ylabel('Velocity [m/s]','fontsize',fs);
    
end