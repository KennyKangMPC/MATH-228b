clear all
pv = [0,0; 1,0; 1,1; 0,1; 0,0];
hmax = 0.15;
nrefmax = 3;
errors = zeros(1, nrefmax);
 
% compute the reference solution
[p_ref, t, e] = pmesh(pv, hmax, nrefmax);
a_ref = fempoi(p_ref, t, e);

flag = 0;
for nref = 0
    [p, t, e] = pmesh(pv, hmax, nref);
    a = fempoi(p, t, e);
    
    % loop over all the points in the ref soln to find the maximum error
    for i = 1:length(p_ref)
        % loop over all the points in the new soln
        for j = 1:length(p)
            % check if x-coordinates match
            if abs(p_ref(i, 1) - p(j, 1)) <= 1e-10
                % check if y-coordinates match
                if abs(p_ref(i, 2) - p(j, 2)) <= 1e-10
                    flag = flag + 1;
                end
            end
        end
    end
end