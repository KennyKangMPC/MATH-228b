% mginit
clear all

% specify the mesh
pv = [0,0; 1,0; 1,1; 0,1; 0,0];
[p, t, e, data, pmid, ia] = pmesh(pv, 0.75, 1);



% % test out the T matrix acting on a linear function z = x + y
% x = data(ref).p(:, 1);
% y = data(ref).p(:, 2);
% funct = x.*x + y.*y; 
% 
% %tplot(data(ref).p, data(ref).t, funct)
% 
% funct_interp = T * funct;
%tplot(data(ref + 1).p, data(ref + 1).t, funct_interp)

%plot(1:length(funct), funct, 'o-', 1:length(funct_interp), funct_interp, '*-')

% perform the solve to get the true solution u0
%[u_true, A, b] = fempoi(p, t, e);
%tplot(p, t, u_true)

%[u, res] = gauss_seidel(A, b, 0 .* u_true, 100);
%semilogy(0:niter, res)

