
l = 6;              % wave number

% generate the mesh
pv = [0,0; 5,0; 5,1; 0,1; 0,0];
[p, t, e] = pmesh(pv, 0.15, 0);

% find the edges of the mesh
[In, Out, Wall] = waveguide_edges(p, t);

LM = t;                         % location matrix                
num_elem = length(LM(:,1));     % number of elements
num_nodes_per_elem = 3;         % linear triangular elements
num_nodes = length(p(:,1));     % total number of nodes

% form the permutation matrix for assembling the global matrices
[perm] = permutation(num_nodes_per_elem);

% specify the boundary conditions - homogeneous neumann and dirichlet
dirichlet_nodes(1,:) = e;       % specify Dirichlet nodes
dirichlet_nodes(2,:) = zeros .* dirichlet_nodes(1,:);
a_k = dirichlet_nodes(2,:);


% specify the quadrature rule
% the rule i calculated
%wt = [1/3, 1/3, 1/3];
%qp = [1/6,1/6; 1/6+3/2,1/6; 1/6,1/6+3/2];

% the one-point rule from last time
wt = [1/2];
qp = [1/3,1/3];

% 1-D quadrature rule
wt1 = [1, 1];
qp1 = [-1/sqrt(3); 1/sqrt(3)];

% two-point rule from last time
%wt = [1/6; 1/6; 1/6];
%qp = [1/6,1/6; 2/3,1/6; 1/6,2/3];

% assemble the elemental k and elemental f
K = zeros(num_nodes); M = zeros(num_nodes);
Bin = zeros(num_nodes); Bout = zeros(num_nodes); F = zeros(num_nodes, 1);

for elem = 1:num_elem
    k = 0;
    f = 0;
    m = 0;

    % compute integrals over the area (entire domain)
    for ll = 1:length(wt)
         [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qp(ll, 1), qp(ll, 2), num_nodes_per_elem, p, LM, elem);
         F_mat = transpose([dx_dxe, dx_deta; dy_dxe, dy_deta]);
         J = det(F_mat);
         D = inv(F_mat) * B;

         % assemble the (elemental) K matrix - correct
         k = k + wt(ll) * transpose(D) * (D) * J;
         
         % assemble the (elemental) m matrix - correct
         m = m + wt(ll) * N * transpose(N) * J;
         
         % assemble the (elemental) forcing vector (b_in)
         f = f + wt(ll) * transpose(N) * J;  
    end
    
    % find the edges of the current element, and sort so that lower number 
    % is first
    edges = [LM(elem,[1,2]); LM(elem,[2,3]); LM(elem,[3,1])];
    edges = sort(edges, 2);
    
    % find if any of the element edges are on the boundaries
    % in In?
    
    for edge = 1:length(edges(:,1))
        % is it on the In boundary?
        for in = 1:length(In(:, 1))
            if edges(edge, 1) == In(in, 1) && edges(edge, 2) == In(in, 2)
                sprintf('Element %i is on In', elem)
            end
        end
        
        % is it on the Out boundary?
        for in = 1:length(Out(:, 1))
            if edges(edge, 1) == Out(in, 1) && edges(edge, 2) == Out(in, 2)
                sprintf('Element %i is on Out', elem)
            end
        end
        
        % is on the Wall boundary?
        for in = 1:length(Wall(:, 1))
            if edges(edge, 1) == Wall(in, 1) && edges(edge, 2) == Wall(in, 2)
                sprintf('Element %i is on Wall', elem)
            end
        end
    end
    
    
    % compute integrals over the boundaries (different quadrature rule)
    % for all boundaries, either xi or eta is constant
    for l = 1:length(wt1)
        % In-boundary (vertical, so either xi or eta is constant)
        [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qp1(l), qp1(l), num_nodes_per_elem, p, LM, elem);
    end
    
    % multiply by 0.5 according to quadrature rule
    k = k .* 0.5;
    m = m .* 0.5;
    f = f .* 0.5;
    
    % place the elemental k matrix into the global K matrix
    for mm = 1:length(perm(:,1))
       i = perm(mm,1);
       j = perm(mm,2);
       K(LM(elem, i), LM(elem, j)) = K(LM(elem, i), LM(elem, j)) + k(i,j);
       M(LM(elem, i), LM(elem, j)) = M(LM(elem, i), LM(elem, j)) + m(i,j);
    end

    % place the elemental f vector into the global F vector
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

%tplot(p, LM, a(1:size(p,1)))