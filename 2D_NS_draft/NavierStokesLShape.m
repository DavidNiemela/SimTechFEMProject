function NavierStokesLShape()

meshSize=0.1;

model = createpde(1);
% standard l-shape geometry is slightly different
geometryFromEdges(model,@lshapeg);

% flip, scale, and translate L-shape geometry to match the required one
model.Geometry = scale(model.Geometry,[-0.5 0.5]);
model.Geometry = translate(model.Geometry,[0.5 0.5]);
generateMesh(model,Hmax=meshSize,GeometricOrder='linear');
[p,e,t] = meshToPet(model.Mesh);

% pdemesh(model.Mesh); % check the mesh

% TODO: for the second problem, can use Channel()
% channel=Channel();
% [p,e,t]=initmesh(channel,'hmax',0.25);

np=size(p,2);
x=p(1,:); y=p(2,:);

% diagonal penalty matrix to enforce zero pressure
out=find(x>0.99); % nodes on outflow
wgts=zeros(np,1); % big weights
wgts(out)=1.e+6;
R=spdiags(wgts,0,np,np);

in = find(y>0.99); % nodes on inflow
bnd=unique([e(1,:) e(2,:)]); % all nodes on boundary
bnd=setdiff(bnd,out); % remove outflow nodes
mask=ones(np,1); % a mask to identify no-slip nodes
mask(bnd)=0; % set mask for no-slip nodes to zero
x=x(in); % x-coordinate of nodes on inflow
y=y(in); % y-coordinate of nodes on inflow

% Dirichlet conditions
g=zeros(np,1); % no-slip values
g(in) = -3 * sin(2*pi * x); % initial profile for y-velocity

[A,~,M]=assema(p,t,1,0,1);
Bx=ConvectionAssembler2D(p,t,ones(np,1),zeros(np,1));
By=ConvectionAssembler2D(p,t,zeros(np,1),ones(np,1));

dt=0.01; % time step
nu=0.1; % viscosity

U=zeros(np,1); % x-velocity
V=zeros(np,1); % y-velocity
P=zeros(np,1); % pressure

for l=1:100
    % enforce no-slip BC
    U=U.*mask;
    V=V.*mask+g;

    % TODO: solve for pressure – is it correct?
    % (for pressure, I simply followed Alg. 29, page 319 of the textbook)
    Pold = P;
    P=(A+R)\-(Bx*U+By*V)/dt;

    % assemble convection matrix
    C=ConvectionAssembler2D(p,t,U,V);
    % intermediate matrices
    K = 0.5*dt*(C + nu*A)./M;
    K1 = eye(np) + K;
    K2 = eye(np) - K;
    % update velocity
    U=K2*U./K1 + 0.5*Bx*(P + Pold);
    V=K2*V./K1 + 0.5*By*(P + Pold);

    pdeplot(p,e,t,'flowdata',[U V]),axis equal,pause(.1)
end

