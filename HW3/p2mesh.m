% first, call pmesh to generate the starting mesh (works for linear
% elements)
clear all

pv = [0,0; 1,0; 1,1; 0,1; 0,0];
hmax = 0.15;
nref = 0.0;
[p,t,e] = pmesh(pv, hmax, nref);g

% plot the starting mesh
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

% scatter(p(:,1), p(:,2), 'o')
% hold on
% scatter(pnew(:,1), pnew(:,2), 'ro')

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
    
    T(i, 2) = row_number(pt1, p2);
    T(i, 4) = row_number(pt2, p2);
    T(i, 6) = row_number(pt3, p2);
    
    % the new boundary nodes are the midpoints between nodes that were
    % previously on the boundary
    
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

% delete entries in e that arise when a triangle has two nodes on a
% boundary, but that triangle spans across the domain such that the
% midpoint is actually inside the domain - this might not work
ecut = [];
for i = 1:length(enew)
    [in, on] = inpolygon(p2(i, 1), p2(i, 2), p(:,1), p(:,2));
    %sprintf('point: %d, in: %i, on: %i', enew(i), in, on)
    if on
        ecut = [ecut; enew(i)];
    end
end

e = [e; ecut];
scatter(p2(e(:,1), 1), p2(e(:,1), 2), 'o')
