% MATH 228b - HW1, question 4
clear all

% number of uniform refinements
nref = 0;

% maximum element size
hmax = 0.2;
hmax = hmax/(2^nref);

% original polygon boundary points
pv = [0,0; 1,0; .5,.5; 1,1; 0,1; 0,0];

% enter sides between which nodes should not connect such that all
% triangles are inside the domain
sides = [2, 3];

% place points on the domain boundaries according to hmax
[new_pts, gmarkers] = initial_mesh(pv, hmax, sides);

% assemble all of the starting node points
pv_orig = pv(1:end-1, :);
pv = [pv_orig; new_pts; 2,2];
plot(pv(:,1), pv(:,2), 'o')

% triangulate the domain
T = delaunayn(pv);

% plot the domain
tplot(pv, T)

% find which triangles are outside the domain, then delete them from T


% find which points are outside the domain, then delete them
[pv] = delete_outside(pv, pv_orig);

% triangulate the domain
T = delaunayn(pv);

% plot the domain
tplot(pv, T)
hold on
for i = 1:2
    plot(pv(T(i,1),1), pv(T(i,1),2), 'o')
    plot(pv(T(i,2),1), pv(T(i,2),2), 'o')
    plot(pv(T(i,3),1), pv(T(i,3),2), 'o')
end


