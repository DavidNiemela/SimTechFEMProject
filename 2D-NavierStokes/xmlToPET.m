function [p,e,t] = xmlToPET(filename, boundaries)
% Takes a xml-meshfile and returns p,e,t Matlab vectors for meshes.
% only for 2d triangle meshes
%input: filename: name of file 'example.xml';
%       boundaries: cellvector containing functionhandles describing points
%                   on outer boundaries
%output: p: vector containing coordinates of points
%        e: vector containing indices of edge-points
%        t: vector defining the nodes of all triangles
meshstruct = readstruct(filename);
mesh = meshstruct.mesh;
t = [vertcat(mesh.cells.triangle.v0Attribute), vertcat(mesh.cells.triangle.v1Attribute), vertcat(mesh.cells.triangle.v2Attribute)]'+1;
x = vertcat(mesh.vertices.vertex.xAttribute)';
y = vertcat(mesh.vertices.vertex.yAttribute)';
p = [x;y];
e = [];
for i = 1:length(boundaries)
    boundaryfunc = boundaries{i};
    b = boundaryfunc(x,y);
    e = [e,b];
end
end

