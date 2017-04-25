function [u, res] = gauss_seidel(A, b, u, niter)
% A - matrix in A * u = b
% b - vector in A * u = b
% u - initial Gauss-Seidel solution guess
% niter - number of Gauss-Seidel iterations

% Gauss-Seidel solution to Au = b
res = zeros(1, niter + 1);
res(1) = max(b - A * u);

for it = 1:niter
    u = u + inv(tril(A)) * (b - A * u);
    res(it + 1) = max(b - A * u);
end

