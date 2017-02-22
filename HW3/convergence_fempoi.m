% script for convergence study of fempoi2
clear all
pv = [0,0; 1,0; 1,1; 0,1; 0,0];
hmax = 0.3;
nrefmax = 4;

errors = zeros(1, nrefmax);
tol = 1e-10;

% create initial (very fine) mesh 
[p_ref, t, e] = pmesh(pv, hmax, nrefmax);
[p2, t2, e2] = p2mesh(p_ref, t);

% compute the reference solution
a_ref = fempoi2(p2, t2, e2);

for nref = 0:(nrefmax - 1)
    [p, t, e] = pmesh(pv, hmax, nref);
    [p, t, e] = p2mesh(p, t);
    a = fempoi2(p, t, e);

    k = 1;
    error = [];

    % loop over all the points in the ref soln to find the maximum error
    for i = 1:length(p_ref)
        % loop over all the points in the new soln
        for j = 1:length(p)
            % check if x-coordinates match
            if abs(p_ref(i, 1) - p(j, 1)) <= tol
                % check if y-coordinates match
                if abs(p_ref(i, 2) - p(j, 2)) <= tol
                    % compute the error at the shared node
                    error(k) = abs(a_ref(i) - a(j));
                    k = k + 1;
                end
            end
        end
    end

    % determine the max-norm (maximum element in error())
    errors(nref + 1) = max(error);
end

loglog(1./hmax .^ (0:3), errors)
for i = 1:(length(errors) - 1)
    rate(i) = log2(errors(i)) - log2(errors(i+1));
end