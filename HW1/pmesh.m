% MATH 228b - HW1, question 4
clear all

% number of uniform refinements
nref = 0;

% maximum element size
hmax = 0.2;
hmax = hmax/(2^nref);

% original polygon boundary points
pv = [0,0; 1,0; .5,.5; 1,1; 0,1; 0,0];

% place points on the domain boundaries according to hmax
[new_pts] = initial_mesh(pv, hmax);

% assemble all of the starting node points
pv_orig = pv(1:end-1, :);
pv = [pv_orig; new_pts; 2,2];
plot(pv(:,1), pv(:,2), 'o')

max = hmax^2 / 2;
area = max + 1;
for r = 1:10
    % find which points are outside the domain, then delete them
    [pv] = delete_outside(pv, pv_orig);
    
    % triangulate the domain
    T = delaunayn(pv);

    % plot the domain
    tplot(pv, T)

    % find which triangles are outside the domain, then delete them from T
    [T] = delete_outside_triangles(T, pv, pv_orig);
    
    % find the circumcenter of the largest triangle
    [pt, max_area] = circumcenter(T, pv);

    % add the circumcenter to the list of points
    pv = [pv; pt];
end





