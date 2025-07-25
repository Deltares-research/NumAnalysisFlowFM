% Gegevens Velema
H     =      80; 
L     =  1800e3;
B     =   550e3;
phi   =      52;
Tmin  =     745;
A0    =       1;

% Basisparameters
T     = Tmin*60;
omega =  2*pi/T;
g     =    9.81; 
aarde = 2*pi/24/3600;
c     = sqrt(g.*H);
f     = 2*aarde*sind(phi);

% Forceringsamplitude
U     =  g/c*A0;

% Rand
cel   =    25e3;
y     = [0:cel:B];

% Inkomende en uitgaande golf
Z1    = c/g*U * exp(-B*f/c) * exp( f/c.*y);
Z2    = c/g*U *     1       * exp(-f/c.*y);

% Fases
phi1  = 0;
phi2  = 45.12;

% Schrijf cmp-files
for i=1:length(y);
    
    % Rooster 1: het grove Cartesische rooster
    namecmp    = ['openboundary1_',num2str(i,'%0.4d'),'.cmp'];
    namecmp    = fopen(namecmp,'wt');
    fprintf(namecmp,['* COLUMNN=3','\n']);
    fprintf(namecmp,['* COLUMN1=Period (min) or Astronomical Componentname','\n']);
    fprintf(namecmp,['* COLUMN2=Amplitude (ISO)','\n']);
    fprintf(namecmp,['* COLUMN3=Phase (deg)','\n']);
    fprintf(namecmp,[num2str(Tmin ,'%7.7f'),'     ',  ...
                     num2str(Z1(i),'%7.7f'),'     ', ...
                     num2str(phi1 ,'%7.7f'),    '\n']);
    fprintf(namecmp,[num2str(Tmin ,'%7.7f'),'     ',  ...
                     num2str(Z2(i),'%7.7f'),'     ', ...
                     num2str(phi2 ,'%7.7f'),    '\n']);
    fclose(namecmp);    
    
    % Rooster 2: het fijne Cartesische rooster
    namecmp    = ['openboundary2_',num2str(i,'%0.4d'),'.cmp'];
    namecmp    = fopen(namecmp,'wt');
    fprintf(namecmp,['* COLUMNN=3','\n']);
    fprintf(namecmp,['* COLUMN1=Period (min) or Astronomical Componentname','\n']);
    fprintf(namecmp,['* COLUMN2=Amplitude (ISO)','\n']);
    fprintf(namecmp,['* COLUMN3=Phase (deg)','\n']);
    fprintf(namecmp,[num2str(Tmin ,'%7.7f'),'     ',  ...
                     num2str(Z1(i),'%7.7f'),'     ', ...
                     num2str(phi1 ,'%7.7f'),    '\n']);
    fprintf(namecmp,[num2str(Tmin ,'%7.7f'),'     ',  ...
                     num2str(Z2(i),'%7.7f'),'     ', ...
                     num2str(phi2 ,'%7.7f'),    '\n']);
    fclose(namecmp);
    
    % Rooster 3: het fijne driehoeksrooster
    namecmp    = ['openboundary3_',num2str(i,'%0.4d'),'.cmp'];
    namecmp    = fopen(namecmp,'wt');
    fprintf(namecmp,['* COLUMNN=3','\n']);
    fprintf(namecmp,['* COLUMN1=Period (min) or Astronomical Componentname','\n']);
    fprintf(namecmp,['* COLUMN2=Amplitude (ISO)','\n']);
    fprintf(namecmp,['* COLUMN3=Phase (deg)','\n']);
    fprintf(namecmp,[num2str(Tmin ,'%7.7f'),'     ',  ...
                     num2str(Z1(i),'%7.7f'),'     ', ...
                     num2str(phi1 ,'%7.7f'),    '\n']);
    fprintf(namecmp,[num2str(Tmin ,'%7.7f'),'     ',  ...
                     num2str(Z2(i),'%7.7f'),'     ', ...
                     num2str(phi2 ,'%7.7f'),    '\n']);
    fclose(namecmp);
end
    