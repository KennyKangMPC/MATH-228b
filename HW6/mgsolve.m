%clear all
function [u, res] = mgsolve(data, vdown, vup, tol, nref)
% data  - structured array containing problem information
% vdown - pre-smoothing G-S iterations
% vup   - post-smoothing G-S iterations
% tol   - res convergence tolerance

% specify the mesh
%pv   = [0,0; 1,0; 1,1; 0,1; 0,0];
%hmax = 0.75;
%nref = 4;

%tol = 1e-10;
%vdown = 100;
%vup = 100;

% precompute all necessary matrices
%data = mginit(pv, hmax, nref);

%nref = fieldCount(data);

% compute the backslash solution for comparison
soln = data(nref + 1).A \ data(nref + 1).b;

res       = [];
loop_iter = 1;
u         = 0 .* data(nref + 1).b;

% compute the first residual
res(loop_iter) = max(abs(data(nref + 1).b - data(nref + 1).A * u));

while res(loop_iter) > tol
    i = nref;
    
    % perform G-S pre-processing to obtain closer solution
    [u] = gauss_seidel(data(i + 1).A, data(i + 1).b, u, vdown);
    
    % begin by solving the actual equation (A * u = b) on the finest mesh
    % compute a residual, then restrict
    r = data(i + 1).b - data(i + 1).A * u;
    r = data(i).R * r;

    % continually solve G-S and then coarsen using recursion. We want to output
    % the error vector of the coarsest mesh using recursion.
    i = i - 1;
    while i > 0
        [e] = gauss_seidel(data(i + 1).A, r, 0 .* r, vdown);

        % compute a residual, then restrict
        r = r - data(i + 1).A * e;
        r = data(i).R * r;
        i = i - 1;
    end

    e = data(i + 1).A \ r;

    % now that we are at the coarsest mesh, we need to interpolate upwards the 
    % error so that we can add it to the solution iterate
    for i = 1:nref
        e = data(i).T * e;
    end

    u = u + e;
    
    % perform G-S post-processing to obtain closer solution
    [u] = gauss_seidel(data(i + 1).A, data(i + 1).b, u, vup);
    
    loop_iter = loop_iter + 1;
    
    % compute the residual for this newest iterate
    res(loop_iter) = max(abs(data(nref + 1).b - data(nref + 1).A * u));
    
end


% plot difference from backslash solution
%tplot(data(nref + 1).p, data(nref + 1).t, abs(u - soln))


end
