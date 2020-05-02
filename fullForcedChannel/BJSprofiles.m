%% BJSprofiles
%% Profiles for channel with BJS boundary conditions

% slip factor
phi = 0.02;

% channle slope
alpha = -0.0005;

% Reynolds
Re = 1e-05;

% Epsilon
eps = 1;

% Nodes in x and z
nx = 2000;
nz = 20;

% Samplig x
xs = 50;

M = csvread('data.csv',1,0);

%% Calculations

% Compute stuff
hend = 1 + alpha*100;
h_in = linspace(0,1,nz);
h_end = linspace(0,hend,nz);

% Height at the sampling point
h = 1+alpha/eps*xs;

% Z coordinate vector at xs

z = linspace(0,h,100);

% Pressure gradient and curvature at xs
dpdx= -3.0/(3.0*phi*h  + h^2);
dpdx2 = (9.0*phi + 6.0*h)*alpha/eps/(3.0*phi*h+h^2)^2;

% Permeation velocity at xs
V = (h + phi)*alpha*z*dpdx - (z.^3/6.0 - h^2*z/2.0 - phi*h*z)*dpdx2;

plot(V,z,'-r')
hold on;
plot(M(:,3),M(:,end),'-b')
    