function [a_orig, K, F] = fempoi2(p, t, e)
%  p = coordinates in mesh
%  t = triangulation
%  e = Dirichlet boundary nodes
%  --- To use this function, you must have already run pmesh()

LM = t;
t2 = t;
p2 = p;
w = 1/3; %quadrature weight
force = 1; %given
coordinates = p;
num_elem = length(LM(:,1)); 

dirichlet_nodes(1,:) = e;       % specify Dirichlet nodes
num_nodes_per_elem = 6;         % quadratic triangular elements

% form the permutation matrix for assembling the global matrices
[perm] = permutation(num_nodes_per_elem);

% specify the boundary conditions - homogeneous neumann and dirichlet
dirichlet_nodes(2,:) = zeros .* dirichlet_nodes(1,:);
a_k = dirichlet_nodes(2,:);
num_nodes = length(p(:,1));

% two-point rule from last time (already multiplied by area)
wt = 0.5.*[1/3; 1/3; 1/3];
qp = [1/6,1/6; 2/3,1/6; 1/6,2/3];

% assemble the elemental k and elemental f
K = zeros(num_nodes);
F = zeros(num_nodes, 1);

for elem = 1:num_elem
    
%     %physical positions of nodes in current element
%     x1 = p2(t2(elem,1),1);
%     x2 = p2(t2(elem,2),1);
%     x3 = p2(t2(elem,3),1);
%     x4 = p2(t2(elem,4),1);
%     x5 = p2(t2(elem,5),1);
%     x6 = p2(t2(elem,6),1);
%     y1 = p2(t2(elem,1),2);
%     y2 = p2(t2(elem,2),2);
%     y3 = p2(t2(elem,3),2);
%     y4 = p2(t2(elem,4),2);
%     y5 = p2(t2(elem,5),2);
%     y6 = p2(t2(elem,6),2);
%     
%     %quadrature points for current element
%     xg1 = (4*x1 + x3 + x5)/6;
%     yg1 = (4*y1 + y3 + y5)/6;
%     xg2 = (x1 + 4*x3 + x5)/6;
%     yg2 = (y1 + 4*y3 + y5)/6;
%     xg3 = (x1 + x3 + 4*x5)/6;
%     yg3 = (y1 + y3 + 4*y5)/6;
%     xg = [xg1 xg2 xg3];
%     yg = [yg1 yg2 yg3];
%     
%     %find coefficients for shape functions
%     X = [x1^2  x2^2  x3^2  x4^2  x5^2  x6^2;
%          y1^2  y2^2  y3^2  y4^2  y5^2  y6^2;
%          x1    x2    x3    x4    x5    x6;
%          y1    y2    y3    y4    y5    y6;
%          x1*y1 x2*y2 x3*y3 x4*y4 x5*y5 x6*y6;
%          1     1     1     1     1     1]; %position matrix
%     C = inv(X); %coefficient matrix
%     
%     %local element area
%     area = polyarea([x1 x3 x5], [y1 y3 y5]);
%     
%     %assemble local stiffness matrix and force vector
%     a = sparse(6,6); %local stiffness matrix
%     bl = sparse(6,1); %local force vector
%     j = 1;
%     while j < 6 + 1
%         %force vector
%         o = 1;
%         while o < 3 + 1 %num quad points + 1
%             bl(j,1) = bl(j,1) + area*w*force*(C(j,1)*xg(o)^2+C(j,2)*yg(o)^2+C(j,3)*xg(o)+C(j,4)*yg(o)+C(j,5)*xg(o)*yg(o)+C(j,6));
%             o = o + 1;
%         end
%             
%         l = 1;
%         while l < 6 + 1
%             %stiffness matrix
%             o = 1;
%             while o < 3 + 1 %num quad points + 1
%                 a(j,l) = a(j,l) + area*w*([2*C(j,1)*xg(o)+C(j,3)+C(j,5)*yg(o), 2*C(j,2)*yg(o)+C(j,4)+C(j,5)*xg(o)]*[2*C(l,1)*xg(o)+C(l,3)+C(l,5)*yg(o); 2*C(l,2)*yg(o)+C(l,4)+C(l,5)*xg(o)]);
%                 o = o + 1;
%             end
%             l = l + 1;
%         end
%         j = j + 1;
%     end
%     
%     k = a;
%     f = bl;
    
    k = 0;
    f = 0;

    for ll = 1:length(wt)
         [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qp(ll, 1), qp(ll, 2), num_nodes_per_elem, coordinates, LM, elem);
         F_mat = transpose([dx_dxe, dx_deta; dy_dxe, dy_deta]);
         J = det(F_mat);

         % assemble the (elemental) forcing vector
         f = f + wt(ll) * transpose(N) * J;

         % assemble the (elemental) stiffness matrix
         k = k + wt(ll) * transpose(inv(F_mat) * B) * inv(F_mat) * B * J;      
    end
    
   if elem == 1
      full(k)
      full(f)
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

% apply Dirichlet conditions
for i = 1:length(e)
    K(e(i),:) = 0;
    K(e(i),e(i)) = 1;
    F(e(i)) = 0;
end

a_orig = K\F;


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

tplot(p, LM, a(1:size(p,1)))

end