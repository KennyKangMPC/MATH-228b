function [a] = fempoi(p, t, e)
%  p = coordinates in mesh
%  t = triangulation
%  e = Dirichlet boundary nodes
%  --- To use this function, you must have already run pmesh()

LM = t;
num_elem = length(LM(:,1)); 
dirichlet_nodes(1,:) = e;
num_nodes_per_elem = 3;         % linear triangular elements

% form the permutation matrix for assembling the global matrices
[perm] = permutation(num_nodes_per_elem);

% specify the boundary conditions - homogeneous neumann and dirichlet
dirichlet_nodes(2,:) = zeros .* dirichlet_nodes(1,:);
a_k = dirichlet_nodes(2,:);
num_nodes = length(p(:,1));

% assemble the elemental k and elemental f
K = zeros(num_nodes);
F = zeros(num_nodes, 1);

for elem = 1:num_elem
    k = 0;
    f = 0;

    X = [p(LM(elem,1),1) p(LM(elem,2),1) p(LM(elem,3),1);
         p(LM(elem,1),2) p(LM(elem,2),2) p(LM(elem,3),2);
         1               1               1               ];
    C = inv(X);
    
    % compute area of triangle based on hard-coded formula
    area = polyarea([p(LM(elem,1),1) p(LM(elem,2),1) p(LM(elem,3),1)], [p(LM(elem,1),2) p(LM(elem,2),2) p(LM(elem,3),2)]);
    
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
       K(LM(elem, i), LM(elem, j)) = K(LM(elem, i), LM(elem, j)) + k(i,j);
    end

    % place the elemental f matrix into the global F matrix
    for i = 1:length(f)
       F(LM(elem, i)) = F((LM(elem, i))) + f(i);
    end
end

% perform static condensation to remove known Dirichlet nodes from solve
[K_uu, K_uk, F_u, F_k] = condensation(K, F, num_nodes, dirichlet_nodes);

% perform the solve
a_u_condensed = K_uu \ (F_u - K_uk * dirichlet_nodes(2,:)');

% expand a_condensed to include the Dirichlet nodes
a = zeros(num_nodes, 1);

a_row = 1;
i = 1;      % index for dirichlet_nodes
j = 1;      % index for expanded row

for a_row = 1:num_nodes
    if (find(dirichlet_nodes(1, :) == a_row))
        a(a_row) = dirichlet_nodes(2,i);
        i = i + 1;
    else
        a(a_row) = a_u_condensed(j);
        j = j + 1;
    end
end

tplot(p, LM, a)

end