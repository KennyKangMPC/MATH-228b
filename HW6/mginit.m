% mginit
clear all

nref = 1;
% specify the mesh
pv = [0,0; 1,0; 1,1; 0,1; 0,0];
[p, t, e, data] = pmesh(pv, 0.75, 1);

% determine the A and b matrices at the finest mesh by calling fempoi
[data(nref + 1).u, data(nref + 1).A, data(nref + 1).b] = ...
    fempoi(data(nref + 1).p, data(nref + 1).t, data(nref + 1).e);

% reduce the finest-mesh A and b 
for i = 1:nref
    data(i).A = data(i).R * data(i + 1).A * data(i).T;
end

%tplot(p, t, u_true)

%[u, res] = gauss_seidel(A, b, 0 .* u_true, 100);
%semilogy(0:niter, res)

