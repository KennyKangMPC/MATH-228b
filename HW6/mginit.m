% mginit

% specify the mesh
pv = [0,0; 2,0; 1.5,1; .5,1; 0,0];
[p, t, e] = pmesh(pv, 0.5, 3);

% perform the solve to get the true solution u0
[u_true, A, b] = fempoi(p, t, e);

[u, res] = gauss_seidel(A, b, 0 .* u_true, niter);
semilogy(0:niter, res)