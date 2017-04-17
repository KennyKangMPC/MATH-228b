
%function [u, res] = mgsolve(data, vdown, vup, tol, nref)
% data  - structured array containing problem information
% vdown - pre-smoothing G-S iterations
% vup   - post-smoothing G-S iterations
% tol   - res convergence tolerance

clear all
pv    = [0,0; 1,0; 1,1; 0,1; 0,0];
hmax  = 0.75;
nref  = 4;
tol   = 1e-10;
vdown = 10;
vup   = 10;
data  = mginit(pv, hmax, nref);

% compute the backslash solution for comparison
soln = data(nref + 1).A \ data(nref + 1).b;

res            = [];
loop_iter      = 1;
u              = 0 .* data(nref + 1).b;
res(loop_iter) = max(abs(data(nref + 1).b - data(nref + 1).A * u));

while res(loop_iter) > tol
    i = nref;
    
    % perform G-S on the original equation (A * u = b)
    [u] = gauss_seidel(data(i + 1).A, data(i + 1).b, u, vdown);
    
    % compute a residual
    r = data(i + 1).b - data(i + 1).A * u;
    
    % restrict to a coarser mesh
    r = data(i).R * r;
    
    % restrict the solution down in prep for upwards v-cycle
    u = data(i).R * u;

    % continually solve G-S and then coarsen using recursion. We want to output
    % the error vector of the coarsest mesh using recursion.
    i = i - 1;
    while i > 0
        % smoothing (initial guess is zero)
        [e] = gauss_seidel(data(i + 1).A, r, 0 .* r, vdown);

        % compute a residual
        r = r - data(i + 1).A * e;
        
        % restrict the residual
        r = data(i).R * r;
        
        % project the solution down in prep for the upwards v-cycle
        u = data(i).R * u;
        
        i = i - 1;
    end

    % analytically solve on the coarsest mesh
    e = data(i + 1).A \ r;

    % now that we are at the coarsest mesh, we need to interpolate upwards the 
    % error so that we can add it to the solution iterate
    for i = 1:nref
        % project back up
        e = data(i).T * e;
        
        % adjust the solution on this mesh
        u = data(i).T * u + e;
        
        % perform G-S iterations to improve projected-up answer
        [u] = gauss_seidel(data(i + 1).A, data(i + 1).b, u, vup);
    end
    
    % perform G-S post-processing to obtain closer solution
    [u] = gauss_seidel(data(i + 1).A, data(i + 1).b, u, vup);
    
    loop_iter = loop_iter + 1;
    
    % compute the residual for this newest iterate
    res(loop_iter) = max(abs(data(nref + 1).b - data(nref + 1).A * u));
    plot(u)
    drawnow
end


% plot difference from backslash solution
%tplot(data(nref + 1).p, data(nref + 1).t, abs(u - soln))


%end
