function B = load_vector_assembler(x,f,u0, uN)
%
% Returns the assembled load vector b.
% Input is a vector x of node coords.
%
N = length(x) - 1;
B = zeros(N+1, 1);
for i = 1:N
h = x(i+1) - x(i);
n = [i i+1];
B(n) = B(n) + [f(x(i)); f(x(i+1))]*h/2;
end
%boundary conditions
B(1) = u0;
B(end) = uN;
