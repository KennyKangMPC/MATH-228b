
function [sol, res] = mgsolve(data, vdown, vup, tol, nref)
% data  - structured array containing problem information
% vdown - pre-smoothing G-S iterations
% vup   - post-smoothing G-S iterations
% tol   - res convergence tolerance

res            = [];
loop_iter      = 1;
data(nref + 1).u = 0 .* data(nref + 1).b;
data(nref + 1).r = data(nref + 1).b - data(nref + 1).A * data(nref + 1).u;
res(loop_iter) = max(abs(data(nref + 1).r));

while res(loop_iter) > tol
    i = nref;
    
    % perform G-S on the original equation (A * u = b)
    [data(i + 1).u] = gauss_seidel(data(i + 1).A, data(i + 1).b, data(i + 1).u, vdown);
    
    % compute a residual (r = b - A * ~u~)
    data(i + 1).r = data(i + 1).b - data(i + 1).A * data(i + 1).u;
    
    % restrict to a coarser mesh (r_h -> r_2h)
    data(i).r = data(i).R * data(i + 1).r;

    for i = fliplr(1:(nref - 1))
        % smoothing the error equation (initial guess is zero)
        % solves (A * e_2h = r_2h)
        [data(i + 1).err] = gauss_seidel(data(i + 1).A, data(i + 1).r, 0 .* data(i + 1).r, vdown);

        % compute a residual (r_2h = r_2h - A * ~e_2h~)
        data(i + 1).r = data(i + 1).r - data(i + 1).A * data(i + 1).err;
        
        % restrict the residual
        data(i).r = data(i).R * data(i + 1).r;
    end

    % analytically solve on the coarsest mesh
    data(1).err = data(1).A \ data(1).r;

    % interpolate upwards
    for i = 1:(nref - 1)
        % project back up
        data(i + 1).err = data(i + 1).err + data(i).T * data(i).err;
  
        % perform G-S iterations to improve projected-up error
        [data(i + 1).err] = gauss_seidel(data(i + 1).A, data(i + 1).r, data(i + 1).err, vup);
    end
    
    % at the last point, the equation is u_h = u_h + e_h;
    data(nref + 1).u = data(nref + 1).u + data(nref).T * data(nref).err;
    
    % smooth u_h
    [data(nref + 1).u] = gauss_seidel(data(nref + 1).A, data(nref + 1).b, data(nref + 1).u, vup);
    
    % compute the residual for this newest iterate
    loop_iter = loop_iter + 1;
    data(nref + 1).r = data(nref + 1).b - data(nref + 1).A * data(nref + 1).u;
    res(loop_iter) = max(abs(data(nref + 1).r));
end

sol = data(nref + 1).u;
end
