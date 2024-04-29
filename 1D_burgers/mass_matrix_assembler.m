function M = mass_matrix_assembler(x)
%
% Returns the assembled stiffness matrix A.
% Input is a vector x of node coords.

N = length(x) - 1; % number of elements
M = zeros(N+1, N+1); % initialize stiffnes matrix to zero
for i = 1:N % loop over elements
h = x(i+1) - x(i); % element length
n = [i i+1]; % nodes
M(n,n) = M(n,n) + [1/3 1/6; 1/6 1/3].*h; % assemble element stiffness
end

%adjust for BC
M(1,1) = 1;
M(1,2) = 0;
M(N+1,N+1) = 1;
M(N+1,N) = 0;
% A(1,1) = 1.e+6;
% A(N+1,N+1) = 1.e+6;