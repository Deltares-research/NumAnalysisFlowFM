function setriemannbndvalues;

format long;
clc;

% Get locations from polyline file
fid         =  fopen('outletBoundary.pli');
data        =  textscan(fid,'%f%f','headerlines',2);
y           =  data{2};
fclose(fid);

% Solid parameters
numseconds  =  23*60*60 + 56*60 + 4.1;
Omega       =  2*pi/numseconds;
g           =  9.81;

% Input
phi         =  45;
u           =  0.1;
d           =  500;
ynul        =  150000;

% Preprocessing
f           =  2.*Omega.*sind(phi);

% Surface elevation
zetafinal   = -f.*(y-ynul)./g.*u;

% Riemann value
H           =  zetafinal + d;
Ri          =  0.5.*(zetafinal + sqrt(H./g).*(-u));

% Write values to files
for i=1:length(y);
    fid     =  fopen(['outletBoundary_',num2str(i,'%0.4d'),'.cmp'],'wt');
    fprintf(fid,'%s\n','* COLUMNN=3');
    fprintf(fid,'%s\n','* COLUMN1=Period (The period of the Riemann invariant)');
    fprintf(fid,'%s\n','* COLUMN2=Amplitude (The amplitude of the Riemann invariant)');
    fprintf(fid,'%s\n','* COLUMN3=Amplitude (The phase of the Riemann invariant)');
    fprintf(fid,'%s',['0.0 ',num2str(Ri(i),'%2.14f'),' 0.0']);
    fclose(fid);
end