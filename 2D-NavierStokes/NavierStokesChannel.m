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
circle = @(x,y) find((x-0.5).^2 + (y-0.5).^2 < (0.1+0.001).^2);
%right boundary
rightboundary = @(x,y) find(x == 3);
boundaries = {leftboundary, topboundary, bottomboundary,circle,rightboundary};

meshtype = 'medium';
meshfile = sprintf('cyl2d_%s.xml', meshtype);
[p,e,t] = xmlToPET(meshfile, boundaries);

np=size(p,2);
x=p(1,:); y=p(2,:);

% diagonal penalty matrix to enforce zero pressure
out = rightboundary(x,y); % nodes on outflow
wgts=zeros(np,1); % big weights
wgts(out)=1.e+6;
R=spdiags(wgts,0,np,np);

in = leftboundary(x,y); % nodes on inflow
% bnd=unique([e(1,:) e(2,:)]); % all nodes on boundary
bnd = e;
bnd=setdiff(bnd,out); % remove outflow nodes
mask=ones(np,1); % a mask to identify no-slip nodes
mask(bnd) = 0; % set mask for no-slip nodes to zero

xin=x(in); % x-coordinate of nodes on inflow
yin=y(in); % y-coordinate of nodes on inflow

% Pmask = ones(np,1);
% Pmask(out) = 0;

% Dirichlet conditions
g=zeros(np,1); % no-slip values
g(in) = 4*sin(pi * yin); % initial profile for y-velocity

A1 = StiffnessAssembler2D(p,t);   % without viscocity
A2 = VariableStiffnessAssembler2D(p,t);  % with viscocity
M = MassAssembler2D(p,t);
Bx=ConvectionAssembler2D(p,t,ones(np,1),zeros(np,1));
By=ConvectionAssembler2D(p,t,zeros(np,1),ones(np,1));

U=zeros(np,1); % x-velocity
V=zeros(np,1); % y-velocity
U=U.*mask+g;
V=V.*mask;

C = ConvectionAssembler2D(p,t,U,V);
P_rhs = Bx*(C*U) + By*(C*V);
P = (A1+R)\(P_rhs);
% P = P.*Pmask;

dt = 0.001; % time step
Tfull = 1;
T = 0;
n = 0;
diviser = fix((Tfull / dt) / 60);
% circle_edge = circle(x,y);
while T < Tfull
    % assemble convection matrix
    C = ConvectionAssembler2D(p,t,U,V);
    
    % constant viscosity
    nu = 0.1;
    LHS_matrix = M - 0.5*dt*(C + nu*A1);
    RHS_matrix = M + 0.5*dt*(C + nu*A1);
    
    % viscosity based on a cell size
    % LHS_matrix = M - (0.5*dt*(C + A2));
    % RHS_matrix = M + 0.5*dt*(C + A2);

    LHS_U = LHS_matrix*U + Bx*P;
    LHS_V = LHS_matrix*V + By*P;

    U = RHS_matrix\LHS_U;
    V = RHS_matrix\LHS_V;

    U=U.*mask+g;
    V=V.*mask;

    C = ConvectionAssembler2D(p,t,U,V);

    P_rhs = Bx*(C*U) + By*(C*V);
    P = (A1+R)\P_rhs;
    % P = P.*Pmask;

    T = T + dt;
    n = n + 1;

    xlim([0 3])
    ylim([0 1])

    if mod(n,diviser) == 0
        % velocity arrows
        quiver(x',y',U,V)
        hold on
        % scatter(x(circle_edge), y(circle_edge), 'r')
        viscircles([0.5,0.5],0.1,'Color','k');
        rectangle('Position',[0 0 3 1])
        pause(0.1)

        title(sprintf('Velocity | mesh: %s | visc: %0.3f', meshtype, nu), 'FontSize', 10);
        text(2.5, 0.8, sprintf('time: %0.2f', T));
        daspect([1 1 1])
        
        fname = sprintf('plots/Channel/%s/U/%d.png', meshtype, n);
        % saveas(gcf, fname, 'png','Resolution',300);
        exportgraphics(gcf,fname,'Resolution',300)
        hold off
    end
end

