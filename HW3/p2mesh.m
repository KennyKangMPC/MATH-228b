% first, call pmesh to generate the starting mesh (works for linear
% elements)
clear all

pv = [0,0; 1,0; 1,1; 0,1; 0,0];
hmax = 0.15;
nref = 0.0;
%[p,t,e] = pmesh(pv, hmax, nref);
p = [0,0; 1,0; 0.5, 1; 1.5,1; 3,0];
t = [1,2,3; 2,4,3; 2,5,4];
e = [1, 2, 3, 4, 5];

% plot the starting mesh
tplot(p, t)
hold on

tol = 1e-10;

% initialize new triangulation
T = zeros(length(t(:,1)), 6);
pnew = [];

% add all of the midpoints
for i = 1:length(t(:,1))
    A = [p(t(i, 1), 1), p(t(i, 1), 2)];
    B = [p(t(i, 2), 1), p(t(i, 2), 2)];
    C = [p(t(i, 3), 1), p(t(i, 3), 2)];

    pnew = [pnew; (A(1) + B(1))/2, (A(2) + B(2))/2];
    pnew = [pnew; (A(1) + C(1))/2, (A(2) + C(2))/2];
    pnew = [pnew; (C(1) + B(1))/2, (C(2) + B(2))/2];
    
    % create new row in T for the triangle - save original points
    T(i, 1) = t(i, 1);
    T(i, 3) = t(i, 2);
    T(i, 5) = t(i, 3);
end

% p2 contains both the original mesh points and all of the midpoints
% (none duplicated)
pnew = unique(pnew, 'rows');
p2 = [p; pnew];

scatter(p(:,1), p(:,2), 'o')
hold on
scatter(pnew(:,1), pnew(:,2), 'ro')

% put in the new points to the triangulation (have to wait until you call
% unique() to make sure the order is consistent

for i = 1:length(t(:,1))
    A = [p(t(i, 1), 1), p(t(i, 1), 2)];
    B = [p(t(i, 2), 1), p(t(i, 2), 2)];
    C = [p(t(i, 3), 1), p(t(i, 3), 2)];
    
    % three midpoint nodes for this triangle
    pt1 = [(A(1) + B(1))/2, (A(2) + B(2))/2];
    pt2 = [(C(1) + B(1))/2, (C(2) + B(2))/2];
    pt3 = [(A(1) + C(1))/2, (A(2) + C(2))/2];
    
    % find the row number of each point in p
    for j = 1:length(p2(:,1))
        % check pt1
        if (abs(pt1(1) - p2(j, 1)) < tol) && (abs(pt1(2) - p2(j, 2)) < tol)
            T(i, 2) = j;
        end
        
        % check pt2
        if (abs(pt2(1) - p2(j, 1)) < tol) && (abs(pt2(2) - p2(j, 2)) < tol)
            T(i, 4) = j;
        end
        
        % check pt3
        if (abs(pt3(1) - p2(j, 1)) < tol) && (abs(pt3(2) - p2(j, 2)) < tol)
            T(i, 6) = j;
        end
    end
end
