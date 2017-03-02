function [ein, eout, ewall] = waveguide_edges(p, t)
tol = 1e-12;

% find all mesh edges - from boundary_nodes.m

edges = [t(:,[1,2]);
         t(:,[2,3]);
         t(:,[3,1])];
edges = sort(edges, 2);
[C, ia, ic] = unique(edges, 'rows');

ein = [];
eout = [];
ewall = [];

for i = 1:length(C)
    flag = 0;
    
    % find vertical edges at x = 0
    if p(C(i, 1), 1) == 0 || p(C(i, 2), 1) == 0
        if abs(p(C(i, 1), 1) - p(C(i, 2), 1)) >= tol
        else
            ein = [ein; C(i, :)];
            flag = 1;
        end
    end
    
    % find vertical edges at x = 5
    if p(C(i, 1), 1) == 5 || p(C(i, 2), 1) == 5
        if abs(p(C(i, 1), 1) - p(C(i, 2), 1)) >= tol
        else
            eout = [eout; C(i, :)];
            flag = 1;
        end
    end
    
    % find horizontal edges with y = 0
    if p(C(i, 1), 2) == 0 || p(C(i, 2), 2) == 0
        if abs(p(C(i, 1), 2) - p(C(i, 2), 2)) >= tol
        else
            ewall = [ewall; C(i, :)];
            flag = 1;
        end
    end
    
    % find horizontal edges with y = 1
    if p(C(i, 1), 2) == 1 || p(C(i, 2), 2) == 1
        if abs(p(C(i, 1), 2) - p(C(i, 2), 2)) >= tol
        else
            ewall = [ewall; C(i, :)];
            flag = 1;
        end
    end
end

% for i = 1:length(ein(:,1))
%     scatter(p(ein(i, 1), 1), p(ein(i, 1), 2), 'ro')
%     hold on
%     scatter(p(ein(i, 2), 1), p(ein(i, 2), 2), 'ro')
%     hold on
% end

% for i = 1:length(eout(:,1))
%     scatter(p(eout(i, 1), 1), p(eout(i, 1), 2), 'bo')
%     hold on
%     scatter(p(eout(i, 2), 1), p(eout(i, 2), 2), 'bo')
%     hold on
% end

% for i = 1:length(ewall(:,1))
%     scatter(p(ewall(i, 1), 1), p(ewall(i, 1), 2), 'go')
%     hold on
%     scatter(p(ewall(i, 2), 1), p(ewall(i, 2), 2), 'go')
%     hold on
% end

