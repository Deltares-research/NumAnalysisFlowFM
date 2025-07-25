function thacker2d;

clc;
close all;
fclose all;
clear all;

% Print?          1 = yes, 0 = no
printen    = 1; 


% Initialisation? 1 = yes, 0 = no
init       = 0;


%%%   INPUT

% Fixed input
g          =      9.81 ;                          
a          = 102000    ; 
D0         =     10    ;
eta        =      2    ; 
L          = 102000    ;


% Derived input
omega      = sqrt(8.*g.*D0)./a;
T          = 2*2*pi/omega;
A          = ((D0 + eta).^2 - D0.^2)./((D0 + eta).^2 + D0.^2);


%%%   WRITE INITIAL SETUP

if init == 1;
    
    %%%   WRITE BED LEVEL BASED ON NET-FILE
    
    % Read x and y
    x          = nc_varget('radial2dsqr_net.nc','NetNode_x');
    y          = nc_varget('radial2dsqr_net.nc','NetNode_y');
    r          = sqrt(x.^2 + y.^2);
    
    % Create z
    z          = -D0.*(1-r.^2./a.^2);
    
    % Write z
    nc_varput('radial2dsqr_net.nc','NetNode_z',z);
    
    
    
    %%%   WRITE WATER LEVEL BASED ON MAP-FILE
    
    % Read x and y
    x          = nc_varget('dflowfmoutput/radial2d_map.nc','FlowElem_xcc');
    y          = nc_varget('dflowfmoutput/radial2d_map.nc','FlowElem_ycc');
    r          = sqrt(x.^2 + y.^2);
    t          = 0;
    
    % Create h
    h          = D0.*(  (sqrt(1-A.^2))./(1-A.*cos(omega.*t))    -  1.0                    ...
                       -(r.^2./L.^2).*((1-A.^2)./(1-A.*cos(omega.*t)).^2  - 1.0)   ) ;
    
    % Write h
    xw         = reshape(x,[size(x,1).*size(x,2) 1]);
    yw         = reshape(y,[size(x,1).*size(x,2) 1]);
    hw         = reshape(h,[size(x,1).*size(x,2) 1]);
    iniwath    = [xw yw hw];
    iniwath    = sortrows(iniwath,2);
    iniwath    = sortrows(iniwath,1);
    dlmwrite('iniwath.xyz', iniwath, 'delimiter', ' ','precision','%7.18f');
    fclose all;

    
else
    
    % Set name
    simname    = ['radial2d'];
    
    % Figure settings
    fs         = 18;
    lw0        =  0.5;
    lw1        =  3;
    lw2        =  2;
    ms         = 14;

    % Read x and y
    clear x y;
    x          = nc_varget(['dflowfmoutput/',simname,'_map.nc'],'FlowElem_xcc');
    y          = nc_varget(['dflowfmoutput/',simname,'_map.nc'],'FlowElem_ycc');
    s          = nc_varget(['dflowfmoutput/',simname,'_map.nc'],'s1');
    d          = nc_varget(['dflowfmoutput/',simname,'_map.nc'],'waterdepth');
    ux         = nc_varget(['dflowfmoutput/',simname,'_map.nc'],'ucx');
    uy         = nc_varget(['dflowfmoutput/',simname,'_map.nc'],'ucy');
    t          = nc_varget(['dflowfmoutput/',simname,'_map.nc'],'time');
    r          = sqrt(x.^2 + y.^2);
    z          = -D0.*(1-r.^2./a.^2);
    un         = sqrt(ux.^2 + uy.^2);
    s(d==0)    = NaN;
    
    % Compute water level at any t
    for j = 1:length(t);
        
        % Define analytical solution
        h          = D0.*(  (sqrt(1-A.^2))./(1-A.*cos(omega.*t(j)))    -  1.0                    ...
                           -(r.^2./L.^2).*((1-A.^2)./(1-A.*cos(omega.*t(j))).^2  - 1.0)   ) ;
        u          = abs(0.5.*r.*omega.*A.*sin(omega.*t(j))/(1-A.*cos(omega.*t(j))));
        z          = -D0.*(1-r.^2./a.^2);
        fa         = find(h<z);
        h(fa)      = NaN;
        u(fa)      = NaN;
        ra         = r;
        ra(fa)     = NaN;
        
        % Catch computational result
        hsim       = s (j,:)';
        usim       = un(j,:)';
        fd         = find(hsim<=z+0.05);
        hsim(fd)   = NaN;
        fd         = isnan(hsim);
        usim(fd)   = NaN;
        rd         = r;
        rd(fd)     = NaN;
        
        % Points for drawing dry/wet interface
        xai        = find(ra==nanmax(ra(:)));
        xa         = ra(xai(1));
        ya(1)      = h(xai(1));
        ya(2)      = u(xai(1));
        xdi        = find(rd==nanmax(rd(:)));
        xd         = rd(xdi(1));
        yd(1)      = hsim(xdi(1));
        yd(2)      = usim(xdi(1));
        
        % Make 1D plot
        figure(2);
        plot(r/1e3,z   ,'k-','linewidth' ,lw1); hold on;
        plot(r/1e3,h   ,'b-','linewidth' ,lw2);
        plot(r/1e3,hsim,'b.','markersize',ms );
        plot(r/1e3,u   ,'r-','linewidth' ,lw2);
        plot(r/1e3,usim,'r.','markersize',ms ); hold off;
        line([xa/1e3 xa/1e3],[ya(1) ya(2)],'linewidth',lw0,'linestyle','--','color','k');
        line([xd/1e3 xd/1e3],[yd(1) yd(2)],'linewidth',lw0,'linestyle','--','color','k');
        grid on;
        
        % Beautify
        xlim([0 120]);
        ylim([-10.0 6]);
        xlabel('radial distance [km]','fontsize',fs);
        ylabel(['orientation [m] / velocity [m/s]'],'fontsize',fs);
        title(['t = ',num2str(t(j)/T*2,'%5.2f'),' \times T, after start (squares)'],'fontsize',fs);
        set(gca,'fontsize',fs);
        
        % Print
        if printen == 1;
            if j == length(t)/2+1;
                legend('geometry','h  analytical','h  D-Flow FM','u  analytical','u  D-Flow FM','location','southeast');
                print(figure(2),'-dpng','-r300','doc/thacker2dsquaresT1.png');
                close all;
            end
            
            if j == length(t);
                legend('geometry','h  analytical','h  D-Flow FM','u  analytical','u  D-Flow FM','location','southeast');
                print(figure(2),'-dpng','-r300','doc/thacker2dsquaresT2.png');
                close all;
            end
            
        end
    
    end
    
end