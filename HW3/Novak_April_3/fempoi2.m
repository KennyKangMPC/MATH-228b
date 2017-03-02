function [a, K, F] = fempoi2(p, t, e)
%  p = p in mesh
%  t = triangulation
%  e = Dirichlet boundary nodes
%  --- To use this function, you must have already run p2mesh()

num_elem = length(t(:,1)); 
num_nodes_per_elem = 6;         % quadratic triangular elements
num_nodes = length(p(:,1));

% form the permutation matrix for assembling the global matrices
[perm] = permutation(num_nodes_per_elem);

% two-point rule from last time (already multiplied by area)
wt = 0.5.*[1/3; 1/3; 1/3];
qp = [1/6,1/6; 2/3,1/6; 1/6,2/3];

% assemble the elemental k and elemental f
K = zeros(num_nodes);
F = zeros(num_nodes, 1);

for elem = 1:num_elem
    k = 0;
    f = 0;

    for ll = 1:length(wt)
         [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qp(ll, 1), qp(ll, 2), num_nodes_per_elem, p, t, elem);
         F_mat = transpose([dx_dxe, dx_deta; dy_dxe, dy_deta]);
         J = det(F_mat);

         % assemble the (elemental) forcing vector
         f = f + wt(ll) * transpose(N) * J;

         % assemble the (elemental) stiffness matrix
         k = k + wt(ll) * transpose(inv(F_mat) * B) * inv(F_mat) * B * J;      
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

% apply Dirichlet conditions
for i = 1:length(e)
    K(e(i),:) = 0;
    K(e(i),e(i)) = 1;
    F(e(i)) = 0;
end

a = K\F;
end