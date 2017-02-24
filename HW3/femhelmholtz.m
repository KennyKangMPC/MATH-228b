
wave = 6;              % wave number

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

% 1-D quadrature rule
wt1 = [1, 1];
qp1 = [-1/sqrt(3); 1/sqrt(3)];

% diagonal edge quadrature rule
wtedge = [1, 1];
qpedge = [-1/sqrt(3), 1--1/sqrt(3); 1/sqrt(3), 1-1/sqrt(3)];

% two-point rule from last time
wt = [1/6; 1/6; 1/6];
qp = [1/6,1/6; 2/3,1/6; 1/6,2/3];

% assemble the elemental k and elemental f
K = zeros(num_nodes); M = zeros(num_nodes);
Bin = zeros(num_nodes); Bout = zeros(num_nodes); F = zeros(num_nodes, 1);

for elem = 1:num_elem
    k = 0; m = 0; bout = 0; bin = 0; bright = 0;

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
    end
    
    % find the edges of the current element
    edges = [LM(elem,[1,2]); LM(elem,[2,3]); LM(elem,[3,1])];
    edges = sort(edges, 2);
    
    in_flag = 0;
    for edge = 1:length(edges(:,1))
        % is it on the In boundary?
        for in = 1:length(In(:, 1))
            if edges(edge, 1) == In(in, 1) && edges(edge, 2) == In(in, 2)
                in_flag = 1;
                % which edge is it - is xe or eta constant?
                xe_edge = [LM(elem, 1), LM(elem, 2)];
                eta_edge = [LM(elem, 3), LM(elem, 1)];
                
                xe_edge = sort(xe_edge, 2);
                eta_edge = sort(eta_edge, 2);
                
                if edges(edge, 1) == xe_edge(1) && edges(edge, 2) == xe_edge(2)
                    for l = 1:length(wt1) % along xe edge (eta = 0)
                        [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qp1(l), 0, num_nodes_per_elem, p, LM, elem);
                        J = abs(p(edges(edge, 1), 2) - p(edges(edge, 2), 2));
                        bin = bin + wt1(l) * N * transpose(N) * J;
                        bright = bright + wt(l) * N * J;
                    end
                elseif edges(edge, 1) == eta_edge(1) && edges(edge, 2) == eta_edge(2)
                    for l = 1:length(wt1) % along eta edge (xe = 0)
                        [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(0, qp1(l), num_nodes_per_elem, p, LM, elem);
                        J = abs(p(edges(edge, 1), 2) - p(edges(edge, 2), 2));
                        bin = bin + wt1(l) * N * transpose(N) * J;
                        bright = bright + wt(l) * N * J;
                    end
                else
                    for l = 1:length(wtedge) % along last edge
                        [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qpedge(l, 1), qpedge(l, 2), num_nodes_per_elem, p, LM, elem);
                        dy = abs(p(edges(edge, 1), 2) - p(edges(edge, 2), 2));
                        dx = abs(p(edges(edge, 1), 1) - p(edges(edge, 2), 1));
                        J = sqrt(dx^2 + dy^2)/sqrt(2);
                        bin = bin + wtedge(l) * N * transpose(N) * J;
                        bright = bright + wt(l) * N * J;
                    end
                end
            end
        end
        
        out_flag = 0;
        % is it on the Out boundary?
        for in = 1:length(Out(:, 1))
            if edges(edge, 1) == Out(in, 1) && edges(edge, 2) == Out(in, 2)
                out_flag = 1;
                % which edge is it - is xe or eta constant?
                xe_edge = [LM(elem, 1), LM(elem, 2)];
                eta_edge = [LM(elem, 3), LM(elem, 1)];
                
                xe_edge = sort(xe_edge, 2);
                eta_edge = sort(eta_edge, 2);
                
                if edges(edge, 1) == xe_edge(1) && edges(edge, 2) == xe_edge(2)
                    for l = 1:length(wt1) % along xe edge (eta = 0)
                        [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qp1(l), 0, num_nodes_per_elem, p, LM, elem);
                        J = abs(p(edges(edge, 1), 2) - p(edges(edge, 2), 2));
                        bout = bout + wt1(l) * N * transpose(N) * J;
                    end
                elseif edges(edge, 1) == eta_edge(1) && edges(edge, 2) == eta_edge(2)
                    for l = 1:length(wt1) % along eta edge (xe = 0)
                        [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(0, qp1(l), num_nodes_per_elem, p, LM, elem);
                        J = abs(p(edges(edge, 1), 2) - p(edges(edge, 2), 2));
                        bout = bout + wt1(l) * N * transpose(N) * J;
                    end
                else
                    for l = 1:length(wtedge) % along last edge
                        [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qpedge(l, 1), qpedge(l, 2), num_nodes_per_elem, p, LM, elem);
                        dy = abs(p(edges(edge, 1), 2) - p(edges(edge, 2), 2));
                        dx = abs(p(edges(edge, 1), 1) - p(edges(edge, 2), 1));
                        J = sqrt(dx^2 + dy^2)/sqrt(2);
                        bout = bout + wtedge(l) * N * transpose(N) * J;
                    end
                end
            end
        end
        
        % is it on the Wall boundary? - do nothing, homogeneous Neumann
    end
    
    % multiply by 0.5 according to quadrature rule
    k = k .* 0.5;
    m = m .* 0.5;
    
    % place the elemental k matrix into the global K matrix
    for mm = 1:length(perm(:,1))
       i = perm(mm,1);
       j = perm(mm,2);
       K(LM(elem, i), LM(elem, j)) = K(LM(elem, i), LM(elem, j)) + k(i,j);
       M(LM(elem, i), LM(elem, j)) = M(LM(elem, i), LM(elem, j)) + m(i,j);
       if (in_flag)
            Bin(LM(elem, i), LM(elem, j)) = Bin(LM(elem, i), LM(elem, j)) + bin(i,j);
       end
       if (out_flag)
           Bout(LM(elem, i), LM(elem, j)) = Bout(LM(elem, i), LM(elem, j)) + bout(i,j);
       end
    end

    % place the elemental f vector into the global F vector
    for i = 1:num_nodes_per_elem
       if (in_flag)
           F(LM(elem, i)) = F((LM(elem, i))) + bright(i);
       end
    end
end

K = K - (wave^2)*M + 1i*wave*(Bin + Bout);
F = F .* 2*1i*wave;

a = K\F;

tplot(p, LM, real(a(1:size(p,1))))