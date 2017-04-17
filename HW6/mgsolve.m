clear all

% specify the mesh
pv   = [0,0; 1,0; 1,1; 0,1; 0,0];
hmax = 0.75;
nref = 4;

niter = 100;
% precompute all necessary matrices
data = mginit(pv, hmax, nref);
% data(i).T promotes from i to i + 1
% data(i).R reduces from i + 1 to i



% begin by solving the actual equation (A * u = b) on the finest mesh
i = nref;
[u] = gauss_seidel(data(i + 1).A, data(i + 1).b, 0 .* data(i + 1).b, niter);

% compute a residual, then restrict
r = data(i + 1).b - data(i + 1).A * u;
r = data(i).R * r;


% continually solve G-S and then coarsen using recursion. We want to output
% the error vector of the coarsest mesh using recursion.
i = i - 1;
%[e] = restrict(data, r, niter, i);

% function that performs some G-S iterations to solve A * u = rhs
% and then computes a residual r = rhs - A * u, and then restricts it

while i > 0
    [e] = gauss_seidel(data(i + 1).A, r, 0 .* r, niter);

    % compute a residual, then restrict
    r = r - data(i + 1).A * e;
    r = data(i).R * r;
    i = i - 1;
end

ge = data(i + 1).A \ r;

% if i == 1
%     e = data(i).A \ r;
%     return
% else
%     [e] = gauss_seidel(data(i + 1).A, r, 0 .* r, niter);
% 
%     % compute a residual, then restrict
%     r = r - data(i + 1).A * e;
%     r = data(i).R * r;
%     i = i - 1;
% end











% solve the error equation with this residual
%[r] = restrict(data(i + 1).A, r, 100, i + 1, data(i).R);





% with this coarser residual, solve the error equation on the coarser mesh,
%e = data(i).A \ r; % but, this is still too difficult to solve, so use G-S
% instead of solving using backslash
% [e, res] = gauss_seidel(data(i).A, r, 0 .* r, 10);
% 
% 
% 
% 
% 
% 
% % then transform this error back to the fine mesh
% e = data(i).T * e;
% 
% % correct the iterative solution using this error
% u = u + e;
% 
% % now, relax again using the updated guess





