close all
clear all

% model = createpde(1);
% standard l-shape geometry is slightly different
% geometryFromEdges(model,@lshapeg);

% flip, scale, and translate L-shape geometry to match the required one
% model.Geometry = scale(model.Geometry,[-0.5 0.5]);
% model.Geometry = translate(model.Geometry,[0.5 0.5]);
% generateMesh(model,Hmax=meshSize,GeometricOrder='linear');


%% READING XML meshes and constructing [p e t]
%edges FOR L SHAPE
%left boundary
leftboundary = @(x,y) find(x==0);
%top boundary
topboundary = @(x,y) find(y == 1);
%bottom boundary
bottomboundary = @(x,y) find(y == 0);
%L boundary
Lboundary = @(x,y) find((y == 0.5 & x>=0.5) | (x == 0.5 & y >= 0.5));
%right boundary
rightboundary = @(x,y) find(x == 1);
boundaries = {leftboundary, topboundary, bottomboundary,Lboundary,rightboundary};

[p,e,t] = xmlToPET('lshape_medium.xml', boundaries);

% [p,e,t] = meshToPet(model.Mesh);

% pdemesh(model.Mesh); % check the mesh

% TODO: for the second problem, can use Channel()
% channel=Channel();
% [p,e,t]=initmesh(channel,'hmax',0.25);

np=size(p,2);
x=p(1,:); y=p(2,:);

% diagonal penalty matrix to enforce zero pressure
out=find(x==1); % nodes on outflow
wgts=zeros(np,1); % big weights
wgts(out)=1.e+6;
R=spdiags(wgts,0,np,np);

in = find(y == 1); % nodes on inflow
% bnd=unique([e(1,:) e(2,:)]); % all nodes on boundary
bnd = e;
bnd=setdiff(bnd,out); % remove outflow nodes
mask=ones(np,1); % a mask to identify no-slip nodes
Pmask = ones(np,1);
mask(bnd)=0; % set mask for no-slip nodes to zero
Pmask(out) = 0;
xin=x(in); % x-coordinate of nodes on inflow
yin=y(in); % y-coordinate of nodes on inflow

% Dirichlet conditions
g=zeros(np,1); % no-slip values
g(in) = -3 * sin(2*pi * xin); % initial profile for y-velocity

A1 = StiffnessAssembler2D(p,t);   % without viscocity
A2 = VariableStiffnessAssembler2D(p,t);  % with viscocity
M = MassAssembler2D(p,t);
% [A,~,M] = assema(p,t,1,0,1);
Bx=ConvectionAssembler2D(p,t,ones(np,1),zeros(np,1));
By=ConvectionAssembler2D(p,t,zeros(np,1),ones(np,1));

dt = 0.0001; % time step


U=zeros(np,1); % x-velocity
V=zeros(np,1); % y-velocity
P=zeros(np,1); % pressure
T = 0;
n = 0;
while T < 1

    U=U.*mask;
    V=V.*mask+g;

    % assemble convection matrix
    C = ConvectionAssembler2D(p,t,U,V);
    
    % constant viscosity
    % nu = 0.001; % viscosity
    % LHS_matrix = M - 0.5*dt*(C + nu*A1);
    % RHS_matrix = M + 0.5*dt*(C + nu*A1);
    
    % % viscosity based on a cell size
    LHS_matrix = M - (0.5*dt*(C + A2));
    RHS_matrix = M + 0.5*dt*(C + A2);

    LHS_U = LHS_matrix*U + Bx*P;    % Shouldn't we multiply Bx and By by M on the left?
    LHS_V = LHS_matrix*V + By*P;

    U = RHS_matrix\LHS_U;
    V = RHS_matrix\LHS_V;
    %
    U=U.*mask;
    V=V.*mask+g;

    C = ConvectionAssembler2D(p,t,U,V);

    P_rhs = Bx*C*U + By*C*V;
    P = (A1+R)\(P_rhs);
    P = P.*Pmask;

    % intermediate matrices
    % K = 0.5*dt*(C + nu*A)./M;
    % K1 = eye(np) + K;
    % K2 = eye(np) - K;
    % update velocity
    % U=K2*U./K1 + 0.5*Bx*(P + Pold);
    % V=K2*V./K1 + 0.5*By*(P + Pold);

    T = T + dt;
    n = n+1;
    % pdeplot(p,e,t,'flowdata',[U V]),axis equal,pause(.1)
    if mod(n,20) == 0
        quiver(x',y',U,V)
        pause(0.1)
    end
end