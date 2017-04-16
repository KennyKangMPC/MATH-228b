function [u, K, F] = fempoi(p, t, e)
%  p = coordinates in mesh
%  t = triangulation
%  e = Dirichlet boundary nodes
%  --- To use this function, you must have already run pmesh()

num_elem = length(t(:,1)); 
dirichlet_nodes(1, :) = e;
num_nodes_per_elem = 3;         % linear triangular elements

% form the permutation matrix for assembling the global matrices
[perm] = permutation(num_nodes_per_elem);

% specify the boundary conditions - homogeneous neumann and dirichlet
dirichlet_nodes(2,:) = zeros .* dirichlet_nodes(1,:);
num_nodes = length(p(:,1));

% assemble the elemental k and elemental f
K = zeros(num_nodes); F = zeros(num_nodes, 1);

for elem = 1:num_elem
    k = 0; f = 0;

    X = [p(t(elem,1),1) p(t(elem,2),1) p(t(elem,3),1);
         p(t(elem,1),2) p(t(elem,2),2) p(t(elem,3),2);
         1              1              1            ];
    C = inv(X);
    
    % compute area of triangle based on hard-coded formula
    area = polyarea([p(t(elem,1),1) p(t(elem,2),1) p(t(elem,3),1)], [p(t(elem,1),2) p(t(elem,2),2) p(t(elem,3),2)]);
    
    for l = 1:3
        f(l, 1) = area / 3;
        for ll = 1:3
            k(l, ll) = (C(ll, 1)*C(l, 1) + C(l, 2)*C(ll, 2)) * area;
        end
    end
    
    % place the elemental k matrix into the global K matrix
    for m = 1:length(perm(:,1))
       i = perm(m,1);
       j = perm(m,2);
       K(t(elem, i), t(elem, j)) = K(t(elem, i), t(elem, j)) + k(i,j);
    end

    % place the elemental f matrix into the global F matrix
    for i = 1:length(f)
       F(t(elem, i)) = F((t(elem, i))) + f(i);
    end
end

% apply Dirichlet boundary conditions
i = 1;
for node = 1:num_nodes
    if (find(dirichlet_nodes(1, :) == node)) % have Dirichlet condition
        for col = 1:num_nodes
            if col == node
                K(node, node) = 1;
            else
                K(node, col) = 0;
            end
        end
        % modify the right-hand side vector - assumes the nodes in 
        % dirichlet_nodes are sorted
        F(node) = dirichlet_nodes(2, i);
        i = i + 1;
    end
end

K = sparse(K);
u = K \ F;
end