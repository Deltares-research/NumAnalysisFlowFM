function [tabel] = computeBedlevelAtFlowlink(mapfile, bedlevtyp);

% clean;
% mapfile           = 'dflowfmoutput\simplechannel_map.nc';
% bedlevtyp         = 3;
% [tabel] = computeBedlevelAtFlowlink('dflowfmoutput\simplechannel_map.nc', 3);

% Input parameters
L        = 10000                  ;
Q        =  2500                  ;
B        =   500                  ;
C        =    65                  ;
deltab   =     1                  ;
ib       = deltab./L              ;
deq      = (Q/B/C/sqrt(ib)).^(2/3);
ueq      = (Q/B/deq)              ;

% Format
format long;

% Tabel van verbinding FlowElem's
FlowLink          = nc_varget(mapfile,'FlowLink'         );

% Coordinaten van flowlink-center
FlowLink_xu       = nc_varget(mapfile,'FlowLink_xu'      );
FlowLink_yu       = nc_varget(mapfile,'FlowLink_yu'      );

% Coordinaten van rondom flowelem's liggende netnodes
FlowElemContour_x = nc_varget(mapfile,'FlowElemContour_x');
FlowElemContour_y = nc_varget(mapfile,'FlowElemContour_y');

% Coordinaten van de cellcenters
FlowElem_xcc      = nc_varget(mapfile,'FlowElem_xcc'     );
FlowElem_ycc      = nc_varget(mapfile,'FlowElem_ycc'     );

% Coordinaten van netnodes
NetNode_x         = nc_varget(mapfile,'NetNode_x'        );
NetNode_y         = nc_varget(mapfile,'NetNode_y'        );
NetNode_z         = nc_varget(mapfile,'NetNode_z'        );

% Waterstand
s1                = nc_varget(mapfile,'s1'               );
s1                = s1(end,:);

% Vind voor elke flowlink de omringende netnodes (slechts 1 van de 2 adjacent cells nodig)
for i=1:size(FlowLink,1);
    
    % Stop indien boundary links behandeld worden
    if FlowLink(i,1) > size(FlowElemContour_x,1) | FlowLink(i,2) > size(FlowElemContour_x,1);
        continue;
    end
    
    % Zet coordinaten van netnodes
    fe1           = FlowLink(i,1);          
    fe2           = FlowLink(i,2);          
    xfl           = FlowElemContour_x(fe1,:);
    yfl           = FlowElemContour_y(fe1,:);
    
    % Kies upwind waterstand (case specifiek!!)
    if FlowElem_xcc(fe1) <= FlowElem_xcc(fe2);
        feupw     = fe1;
    else
        feupw     = fe2;
    end
    
    % Check voor elke coordinaat welke netnode het is
    for j=1:length(xfl);
        netnodes  = [NetNode_x NetNode_y];
        [ii,jj]   = find(netnodes(:,1)==xfl(j) & netnodes(:,2)==yfl(j));
        nnx(j)    = ii;
        nny(j)    = ii;
    end
    nnx(j+1)      = nnx(1);
    nny(j+1)      = nny(1);
    
    % Check duo's van netnodes
    for k=1:length(nnx)-1;
        xn1       = NetNode_x(nnx(k  ));
        xn2       = NetNode_x(nnx(k+1));
        yn1       = NetNode_y(nny(k  ));
        yn2       = NetNode_y(nny(k+1));
        xmean     = 0.5.*(xn1 + xn2);
        ymean     = 0.5.*(yn1 + yn2);
        if (xmean == FlowLink_xu(i) & ymean == FlowLink_yu(i));
            netnode1       = min(nnx(k  ),nnx(k+1));
            netnode2       = max(nnx(k  ),nnx(k+1));
            tabel(i,1)     = i;
            tabel(i,2)     = netnode1;
            tabel(i,3)     = netnode2;
            tabel(i,4)     = fe1;
            tabel(i,5)     = fe2;
            tabel(i,6)     = NetNode_z(netnode1);
            tabel(i,7)     = NetNode_z(netnode2);
            if     bedlevtyp == 3;
                zk(i)      = 0.5.*(NetNode_z(netnode1) + NetNode_z(netnode2));
            elseif bedlevtyp == 4;
                zk(i)      = min  (NetNode_z(netnode1),  NetNode_z(netnode2));
            elseif bedlevtyp == 5;
                zk(i)      = max  (NetNode_z(netnode1),  NetNode_z(netnode2));
            end
            tabel(i,8)     = s1(feupw);
            tabel(i,9)     = zk(i);
            tabel(i,10)    = FlowLink_xu(i);
            tabel(i,11)    = FlowLink_yu(i);
            tabel(i,12)    = s1(feupw) - zk(i);
            continue
        end
    end
    
    % Ruim op
    clear nnx nny;
end