%function [p2, t2, e2] = p2mesh(p, t, e)
clear all
% plot the starting mesh
% tplot(p, t)
% hold on
%pv = [0,0; 1,0; 1,1; 0,1; 0,0];
pv = [0,0; 1,0; 0.5,0.5; 1,1; 0,1; 0,0];
[p, t, e] = pmesh(pv, 0.15, 0.0);

e = boundary_nodes(t);


%p = [0,0; 1,0; 0.5,1; 1.5,1; 3,0];
%t = [1,2,3; 2,4,3; 2,5,4];
%e = [1,2,3,4,5];
tplot(p, t)
hold on

% initialize new triangulation and the new nodes and boundary nodes
T = zeros(length(t(:,1)), 6);
pnew = [];
enew = [];

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
scatter(p2(:,1), p2(:,2), 'bo')
hold on

for i = 1:length(t(:,1))
    A = [p(t(i, 1), 1), p(t(i, 1), 2)];
    B = [p(t(i, 2), 1), p(t(i, 2), 2)];
    C = [p(t(i, 3), 1), p(t(i, 3), 2)];
    
    % three midpoint nodes for this triangle
    pt1 = [(A(1) + B(1))/2, (A(2) + B(2))/2];
    pt2 = [(C(1) + B(1))/2, (C(2) + B(2))/2];
    pt3 = [(A(1) + C(1))/2, (A(2) + C(2))/2];
    
    T(i, 2) = row_number(pt1, p2);
    T(i, 4) = row_number(pt2, p2);
    T(i, 6) = row_number(pt3, p2);
    
    % check if side 1-2 or 1-3 is on the boundary (both nodes in e)
    for j = 1:length(e)
        if e(j) == T(i, 1) % the 1-coordinate is on the boundary
            for k = 1:length(e)
                % either the 2 or 3-coordinate is on boundary
                if e(k) == T(i, 3)
                    %sprintf('1: %i, 2: %i', T(i, 1), T(i, 3))
                    enew = [enew; row_number(pt1, p2)];
                end
                if e(k) == T(i, 5)
                    enew = [enew; row_number(pt3, p2)];
                end
            end
        end 
    end
    
    % check if side 2-3 is on the boundary (both nodes in e)
    for j = 1:length(e)
        if e(j) == T(i, 3) % the 2-coordinate is on the boundary
            for k = 1:length(e)
                if e(k) == T(i, 5) % the 3-coordinate is on the boundary
                    enew = [enew; row_number(pt2, p2)];
                end
            end
        end
    end
end

enew = unique(enew);
t2 = T;

% delete points in e2 that appear more than one time in the triangulation
enewer = [];
for k = 1:length(enew)
    flag = 0;

    for i = 1:length(T(:,1))
        for j = 1:6
            if T(i, j) == enew(k)
                flag = flag + 1;
            end
        end
    end
    
    if flag >= 2
        % point isn't actually on boundary
    elseif flag == 1
        enewer = [enewer; enew(k)];
    end
end

e2 = [e; enewer];

for i = 1:length(e2)
    scatter(p2(e2(i), 1), p2(e2(i), 2), 'ro')
    hold on
end
