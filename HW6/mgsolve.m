clear all

% specify the mesh
pv   = [0,0; 1,0; 1,1; 0,1; 0,0];
hmax = 0.75;
nref = 1;

% precompute all necessary matrices
data = mginit(pv, hmax, nref);
% data(i).T promotes from i to i + 1
% data(i).R reduces from i + 1 to i

% perform first G-S iterations (10 relaxations)
[u, res] = gauss_seidel(data(nref + 1).A, data(nref + 1).b, 0 .* data(nref + 1).b, 10);

% compute a residual
r = data(nref + 1).b - data(nref + 1).A * u;

% restrict
r = data(nref).R * r;

% with this coarser residual, solve the error equation on the coarser mesh
e = data(nref).A \ r;

% transform this error back to the fine mesh
e = data(nref).T * e;

% correct the iterative solution using this error
u = u + e;

% now, relax again using the updated guess

