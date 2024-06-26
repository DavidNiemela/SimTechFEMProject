close all
clear all

%% READING XML meshes and constructing [p e t]
%left boundary
leftboundary = @(x,y) find(x==0);
%top boundary
topboundary = @(x,y) find(y == 1);
%bottom boundary
bottomboundary = @(x,y) find(y == 0);
%L boundary
Lboundary = @(x,y) find((y == 0.5 & x>=0.5) | (x == 0.5 & y >= 0.5));
%right boundary
rightboundary = @(x,y) find((x == 1 & y ~= 0) & (x == 1 & y ~= 0.5));
boundaries = {leftboundary, topboundary, bottomboundary,Lboundary,rightboundary};

meshtype = 'coarse';
meshfile = sprintf('lshape_%s.xml', meshtype);
[p,e,t] = xmlToPET(meshfile, boundaries);

np=size(p,2);
x=p(1,:); y=p(2,:);

% diagonal penalty matrix to enforce zero pressure
out = rightboundary(x,y); % nodes on outflow
wgts = zeros(np,1); % big weights
wgts(out) = 1.e+6;
R=spdiags(wgts,0,np,np);

in = topboundary(x,y); % nodes on inflow
% bnd=unique([e(1,:) e(2,:)]); % all nodes on boundary
bnd = e;
bnd = setdiff(bnd, out); % remove outflow nodes
mask = ones(np,1); % a mask to identify no-slip nodes
mask(bnd) = 0; % set mask for no-slip nodes to zero

xin=x(in); % x-coordinate of nodes on inflow
yin=y(in); % y-coordinate of nodes on inflow

Pmask = ones(np,1);
Pmask(out) = 0;

% Dirichlet conditions
g = zeros(np,1); % no-slip values
g(in) = -3 * sin(2*pi * xin); % initial profile for y-velocity

A1 = StiffnessAssembler2D(p,t);   % without viscocity
A2 = VariableStiffnessAssembler2D(p,t);  % with viscocity
M = MassAssembler2D(p,t);
Bx = ConvectionAssembler2D(p,t,ones(np,1),zeros(np,1));
By = ConvectionAssembler2D(p,t,zeros(np,1),ones(np,1));

U=zeros(np,1); % x-velocity
V=zeros(np,1); % y-velocity
U=U.*mask;
V=V.*mask+g;

C = ConvectionAssembler2D(p,t,U,V);
P_rhs = Bx*C*U + By*C*V;
P = (A1+R)\(P_rhs);
% P = P.*Pmask;

%% Specify viscosity (0 for mesh dependent)
nu = 0.01;

%% Create folders for plots
folder = sprintf('plots/LShape/%s/nu_%0.3f', meshtype, nu);
if ~exist(folder, 'dir')
    mkdir(folder)
end
plots = {'U', 'Uabs', 'P'};
for i = 1:length(plots)
    varFolder = strcat(folder,sprintf('/%s', plots{i}));
    if ~exist(varFolder, 'dir')
        mkdir(varFolder)
    end
end

%% Time steps
dt = 0.01; % time step
Tfull = 1;
T = 0;
n = 0;
diviser = fix((Tfull / dt) / 60);
while T < 1
    % assemble convection matrix
    C = ConvectionAssembler2D(p,t,U,V);
    
    if nu > 0
        % constant viscosity
        LHS_matrix = M - 0.5*dt*(C + nu*A1);
        RHS_matrix = M + 0.5*dt*(C + nu*A1);
    else
        % viscosity based on a cell size
        LHS_matrix = M - (0.5*dt*(C + A2));
        RHS_matrix = M + 0.5*dt*(C + A2);
    end

    LHS_U = LHS_matrix*U + Bx*P;
    LHS_V = LHS_matrix*V + By*P;

    U = RHS_matrix\LHS_U;
    V = RHS_matrix\LHS_V;

    U=U.*mask;
    V=V.*mask+g;

    C = ConvectionAssembler2D(p,t,U,V);

    P_rhs = Bx*(C*U) + By*(C*V);
    P = (A1+R)\P_rhs; 
    % P = P.*Pmask;

    T = T + dt;
    n = n + 1;
    % pdeplot(p,e,t,'flowdata',[U V]),axis equal,pause(.1)
    if mod(n,diviser) == 0
        plotVelocityInLShape(x,y,U,V,meshtype,nu,T);
        plotAbsVelocityInLShape(p,t,U,V,meshtype,nu,T);
        plotPressureInLShape(p,t,P,meshtype,nu,T);
    end
end