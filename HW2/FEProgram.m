clear all

num_nodes_per_elem = 3;         % linear triangular elements

% form the permutation matrix for assembling the global matrices
[permutation] = permutation(num_nodes_per_elem);

% perform meshing using unstructured triangles
% t represents the Location Matrix (LM)
[p, LM, e] = pmesh([0,0; 1,0; 1,1; 0,1; 0,0], 0.15, 0);
num_elem = length(LM(:,1)); 

% specify the boundary conditions - homogeneous neumann and dirichlet
dirichlet_nodes(1,:) = e(p(e,1) < 1e-6 | p(e,2) < 1e-6);
dirichlet_nodes(2,:) = zeros .* dirichlet_nodes(1,:);
a_k = dirichlet_nodes(2,:);
num_nodes = length(p(:,1));

% define the quadrature rule
[wt, qp] = quadrature(num_nodes_per_elem);

% assemble the elemental k and elemental f
K = zeros(num_nodes);
F = zeros(num_nodes, 1);

for elem = 1:num_elem
    k = 0;
    f = 0;

    for ll = 1:length(qp) % eta loop
         for l = 1:length(qp) % xe loop
            [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qp(l), qp(ll), num_nodes_per_elem, p, LM, elem);
            F_mat = transpose([dx_dxe, dx_deta; dy_dxe, dy_deta]);
            J = det(F_mat);

            % assemble the (elemental) forcing vector
            f = f + wt(ll) * wt(l) * transpose(N) * J;

            % assemble the (elemental) stiffness matrix - correct
            k = k + wt(ll) * wt(l) * transpose(inv(F_mat) * B) * inv(F_mat) * B * J;      
         end
    end

    % place the elemental k matrix into the global K matrix
    for m = 1:length(permutation(:,1))
       i = permutation(m,1);
       j = permutation(m,2);
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



