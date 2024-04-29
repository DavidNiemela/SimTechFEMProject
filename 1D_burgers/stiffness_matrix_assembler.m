function A = stiffness_matrix_assembler(x)
%
% Returns the assembled stiffness matrix A.
% Input is a vector x of node coords.

N = length(x) - 1; % number of elements
A = zeros(N+1, N+1); % initialize stiffnes matrix to zero
for i = 1:N % loop over elements
h = x(i+1) - x(i); % element lengthfunc
n = [i i+1]; % nodes
A(n,n) = A(n,n) + [1 -1; -1 1]/h; % assemble element stiffness
end

%adjust for BC
A(1,1) = 1;
A(1,2) = 0;
A(N+1,N+1) = 1;
A(N+1,N) = 0;
% A(1,1) = 1.e+6;
% A(N+1,N+1) = 1.e+6;