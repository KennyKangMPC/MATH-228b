clear all

% specify the mesh
pv   = [0,0; 1,0; 1,1; 0,1; 0,0];
hmax = 0.75;
nref = 4;

% precompute all necessary matrices
data = mginit(pv, hmax, nref);
% data(i).T promotes from i to i + 1
% data(i).R reduces from i + 1 to i


g


% obtain the residual needed for the first error equation
i = nref;
[r] = restrict(data(i + 1).A, data(i + 1).b, 100, i + 1, data(i).R);

% solve the error equation with this residual
%[r] = restrict(data(i).A, r, 100, i, data(i - 1).R);





% with this coarser residual, solve the error equation on the coarser mesh,
%e = data(i).A \ r; % but, this is still too difficult to solve, so use G-S
% instead of solving using backslash
[e, res] = gauss_seidel(data(i).A, r, 0 .* r, 10);






% then transform this error back to the fine mesh
e = data(i).T * e;

% correct the iterative solution using this error
u = u + e;

% now, relax again using the updated guess





