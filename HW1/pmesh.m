% MATH 228b - HW1, question 4
clear all

% number of uniform refinements
nref = 0;

% maximum element size
hmax = 0.2;
hmax = hmax/(2^nref);

% original polygon
pv = [0,0; 1,0; .5,.5; 1,1; 0,1; 0,0];

% place points on the domain boundaries according to hmax
[new_pts] = initial_mesh(pv, hmax);

% assemble all of the starting node points
pv_orig = pv(1:end-1, :);
pv = [pv_orig; new_pts; 2,2];
plot(pv(:,1), pv(:,2), 'o')

% triangulate the domain
T = delaunayn(pv);

% plot the domain
tplot(pv, T)

j = 1;
% find which points are outside the domain, then delete them
for i = 1:size(pv(:,1))
    in = inpolygon(pv(i,1), pv(i,2), pv_orig(:,1), pv_orig(:,2));
    
    if in == 0
        % generic point outside the domain
    else
        % inside the domain
        pv_inside(j,:) = pv(i, :);
        j = j + 1;
    end
end

pv = pv_inside;

% triangulate the domain
T = delaunayn(pv);

% plot the domain
tplot(pv, T)