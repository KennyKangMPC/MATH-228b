% MATH 228b - HW1, question 4
clear all

% number of uniform refinements
nref = 2;

% maximum element size
hmax = 0.4;

% original polygon boundary points
pv = [0,0; 1,0; .5,.5; 1,1; 0,1; 0,0];
%pv = [0,0; 1,0; 2,1; 2,2; 1,3; 0,3; -1,2; -1,1; 0,0];
pv = [2,0; 3,1; 2,2; 0,2; 2,0];

% place points on the domain boundaries according to hmax
[new_pts] = initial_mesh(pv, hmax);

% assemble all of the starting node points
pv_orig = pv(1:end-1, :);
pv = [pv_orig; new_pts; 2,2];

max = hmax^2 / 2;
max_area = max + 1; % dummy initial value
while max_area > max
    % find which points are outside the domain, then delete them
    pv = unique(pv, 'rows');
    [pv] = delete_outside(pv, pv_orig);
    
    % triangulate the domain
    T = delaunayn(pv);

    % find which triangles are outside the domain, then delete them from T
    [T] = delete_outside_triangles(T, pv, pv_orig);
    tplot(pv, T)
    
    % compute triangle areas using Heron's formula
    for i = 1:length(T(:,1))
        A = [pv(T(i, 1), 1), pv(T(i, 1), 2)];
        B = [pv(T(i, 2), 1), pv(T(i, 2), 2)];
        C = [pv(T(i, 3), 1), pv(T(i, 3), 2)];
        area(i) = (A(1) * (B(2) - C(2)) + B(1) * (C(2) - A(2)) + C(1) * (A(2) - B(2))) / 2;
        if i == 1
            max_area = area(1);
        end
        
        if area(i) > max_area
            max_area = area(i);
            refine = i;
        end 
    end

    A = [pv(T(refine, 1), 1), pv(T(refine, 1), 2)];
    B = [pv(T(refine, 2), 1), pv(T(refine, 2), 2)];
    C = [pv(T(refine, 3), 1), pv(T(refine, 3), 2)];
    
    [pt] = circumcenter(A, B, C, pv);

    % add the circumcenter to the list of points
    pv = [pv; pt];
end

% remove last triangulation (not needed since we just broke from loop)
pv = pv(1:(end-1), :);

tplot(pv, T)
saveas(gcf, 'refine-B0', 'png')

% perform any uniform refinements
for uf = 1:nref
    sprintf('Performing uniform refinement %i', uf)

    for i = 1:length(T(:,1))
         A = [pv(T(i, 1), 1), pv(T(i, 1), 2)];
         B = [pv(T(i, 2), 1), pv(T(i, 2), 2)];
         C = [pv(T(i, 3), 1), pv(T(i, 3), 2)];
        
        % add the three midpoints
        pv = [pv; (A(1) + B(1))/2, (A(2) + B(2))/2];
        pv = [pv; (A(1) + C(1))/2, (A(2) + C(2))/2];
        pv = [pv; (C(1) + B(1))/2, (C(2) + B(2))/2];
    end
    
    pv = unique(pv, 'rows');
    T = delaunayn(pv);
    
    % find which triangles are outside the domain, then delete them from T
    [T] = delete_outside_triangles(T, pv, pv_orig);
    
    tplot(pv, T)
    saveas(gcf, sprintf('refine-B%i',uf), 'png');

end


