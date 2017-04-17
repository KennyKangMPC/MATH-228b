% mginit
clear all

nref = 1;
% specify the mesh
pv = [0,0; 1,0; 1,1; 0,1; 0,0];
[p, t, e, data] = pmesh(pv, 0.75, 1);

% determine the A and b matrices by calling fempoi
for i = 1:(nref + 1)
    [data(i).u, data(i).A, data(i).b] = fempoi(data(i).p, data(i).t, data(i).e);
end
% perform the solve to get the true solution u0

%tplot(p, t, u_true)

%[u, res] = gauss_seidel(A, b, 0 .* u_true, 100);
%semilogy(0:niter, res)

