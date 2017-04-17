function [r] = restrict(data, rhs, niter, i)
% function that performs some G-S iterations to solve A * u = rhs
% and then computes a residual r = rhs - A * u, and then restricts it

% If the mesh is already at its finest point, then break from recursion

if i == 1
    r = data(i).A \ rhs;
    return
else
    [e] = gauss_seidel(data(i + 1).A, rhs, 0 .* rhs, niter);

    % compute a residual, then restrict
    r = rhs - data(i + 1).A * e;
    r = data(i).R * r;
    rhs = r; % rhs becomes this computed residual
    i = i - 1;
end

