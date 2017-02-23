function [In, Out, Wall] = waveguide_edges(p, t)
tol = 1e-12;

% find all mesh edges - from boundary_nodes.m

edges = [t(:,[1,2]);
         t(:,[2,3]);
         t(:,[3,1])];
edges = sort(edges, 2);
[C, ia, ic] = unique(edges, 'rows');

In = [];
Out = [];
Wall = [];

for i = 1:length(C)
    flag = 0;
    
    % find vertical edges at x = 0
    if p(C(i, 1), 1) == 0 || p(C(i, 2), 1) == 0
        if abs(p(C(i, 1), 1) - p(C(i, 2), 1)) >= tol
        else
            In = [In; C(i, :)];
            flag = 1;
        end
    end
    
    % find vertical edges at x = 5
    if p(C(i, 1), 1) == 5 || p(C(i, 2), 1) == 5
        if abs(p(C(i, 1), 1) - p(C(i, 2), 1)) >= tol
        else
            Out = [Out; C(i, :)];
            flag = 1;
        end
    end
    
    % find horizontal edges with y = 0
    if p(C(i, 1), 2) == 0 || p(C(i, 2), 2) == 0
        if abs(p(C(i, 1), 2) - p(C(i, 2), 2)) >= tol
        else
            Wall = [Wall; C(i, :)];
            flag = 1;
        end
    end
    
    % find horizontal edges with y = 1
    if p(C(i, 1), 2) == 1 || p(C(i, 2), 2) == 1
        if abs(p(C(i, 1), 2) - p(C(i, 2), 2)) >= tol
        else
            Wall = [Wall; C(i, :)];
            flag = 1;
        end
    end
end

for i = 1:length(In(:,1))
    scatter(p(In(i, 1), 1), p(In(i, 1), 2), 'ro')
    hold on
    scatter(p(In(i, 2), 1), p(In(i, 2), 2), 'ro')
    hold on
end

for i = 1:length(Out(:,1))
    scatter(p(Out(i, 1), 1), p(Out(i, 1), 2), 'bo')
    hold on
    scatter(p(Out(i, 2), 1), p(Out(i, 2), 2), 'bo')
    hold on
end

for i = 1:length(Wall(:,1))
    scatter(p(Wall(i, 1), 1), p(Wall(i, 1), 2), 'go')
    hold on
    scatter(p(Wall(i, 2), 1), p(Wall(i, 2), 2), 'go')
    hold on
end

