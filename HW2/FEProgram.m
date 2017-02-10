clear all

fontsize = 16;                  % fontsize for plots
num_nodes_per_elem = 3;         % linear triangular elements

% form the permutation matrix for assembling the global matrices
[permutation] = permutation(num_nodes_per_elem);



   
%     % for a 2-D mesh
%     [coordinates, LM] = polar_mesh(No, Nr, dt, num_nodes, ri, ro, num_elem);
    
%     % specify the boundary conditions
%     [dirichlet_nodes, neumann_nodes, a_k] = BCnodes(left, right, left_value, right_value, num_nodes, Nr, No);
% 
%     % define the quadrature rule
%     [wt, qp] = quadrature(shape_order);
% 
%     % assemble the elemental k and elemental f
%     K = zeros(num_nodes);
%     F = zeros(num_nodes, 1);
% 
%     for elem = 1:num_elem
%         k = 0;
%         f = 0;
% 
%         for ll = 1:length(qp) % eta loop
%              for l = 1:length(qp) % xe loop
%                      [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qp(l), qp(ll), num_nodes_per_elem, coordinates, LM, elem);
%                      F_mat = transpose([dx_dxe, dx_deta; dy_dxe, dy_deta]);
%                      J = det(F_mat);
%                      r = sqrt(x_xe_eta^2 + y_xe_eta^2);
%                      
%                      % assemble the (elemental) forcing vector
%                      f = f + wt(ll) * wt(l) * (40 * sin(2 * theta) / (r^2)) * transpose(N) * J;
%                      
%                      % assemble the (elemental) stiffness matrix - correct
%                      k = k + wt(ll) * wt(l) * transpose(inv(F_mat) * B) * k_th * inv(F_mat) * B * J;      
%              end
%              
%         end
% 
%          % place the elemental k matrix into the global K matrix
%          for m = 1:length(permutation(:,1))
%             i = permutation(m,1);
%             j = permutation(m,2);
%             K(LM(elem, i), LM(elem, j)) = K(LM(elem, i), LM(elem, j)) + k(i,j);
%          end
% 
%          % place the elemental f matrix into the global F matrix
%          for i = 1:length(f)
%             F(LM(elem, i)) = F((LM(elem, i))) + f(i);
%          end
%     end
% 
% % perform static condensation to remove known Dirichlet nodes from solve
% [K_uu, K_uk, F_u, F_k] = condensation(K, F, num_nodes, dirichlet_nodes);
% 
% % perform the solve
% a_u_condensed = K_uu \ (F_u - K_uk * dirichlet_nodes(2,:)');
% 
% % expand a_condensed to include the Dirichlet nodes
% a = zeros(num_nodes, 1);
% 
% a_row = 1;
% i = 1;      % index for dirichlet_nodes
% j = 1;      % index for expanded row
% 
% for a_row = 1:num_nodes
%     if (find(dirichlet_nodes(1, :) == a_row))
%         a(a_row) = dirichlet_nodes(2,i);
%         i = i + 1;
%     else
%         a(a_row) = a_u_condensed(j);
%         j = j + 1;
%     end
% end
% 
% % assemble the solution in the physical domain
% [mat] = postprocess(num_elem, parent_domain, a, LM, num_nodes_per_elem, shape_order, coordinates, physical_domain);



