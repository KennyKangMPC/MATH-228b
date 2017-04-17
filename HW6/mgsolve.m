clear all

% specify the mesh
pv   = [0,0; 1,0; 1,1; 0,1; 0,0];
hmax = 0.1;
nref = 1;

% precompute all necessary matrices
data = mginit(pv, hmax, nref);