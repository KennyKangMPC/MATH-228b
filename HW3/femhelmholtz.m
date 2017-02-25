function [K, M, Bin, Bout, bin] = femhelmholtz(p, t, ein, eout)
wave = 6;              % wave number
imag = sqrt(-1);
               
num_elem = length(t(:,1));      % number of elements
num_nodes_per_elem = 3;         % linear triangular elements
num_nodes = length(p(:,1));     % total number of nodes

% form the permutation matrix for assembling the global matrices
[perm] = permutation(num_nodes_per_elem);

% 1-D quadrature rule (Simpson's rule)
wt1 = [1/6, 4/6, 1/6];
qp1 = [0, 0.5, 1];

% two-point rule - already multiplied by the area
wt = 0.5 .* [1/3; 1/3; 1/3];
qp = [1/6,1/6; 2/3,1/6; 1/6,2/3];

K = sparse(num_nodes, num_nodes); M = sparse(num_nodes, num_nodes);
Bin = sparse(num_nodes, num_nodes); Bout = sparse(num_nodes, num_nodes); F = zeros(num_nodes, 1);

% change the t for weird elements
[t] = changeLM(num_elem, t, ein, eout);

for elem = 1:num_elem
    k = 0; m = 0; bout = 0; bin = 0; bright = 0;

    for ll = 1:length(wt)
         [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qp(ll, 1), qp(ll, 2), num_nodes_per_elem, p, t, elem);
         F_mat = transpose([dx_dxe, dx_deta; dy_dxe, dy_deta]);
         J = det(F_mat);
         D = inv(F_mat) * B;

         k = k + wt(ll) * transpose(D) * (D) * J;
         m = m + wt(ll) * N * transpose(N) * J; 
    end
    
    % find the edges of the current element
    edges = [t(elem,[1,2]); t(elem,[2,3]); t(elem,[3,1])];
    edges = sort(edges, 2);
    
    xe_edge = [t(elem, 1), t(elem, 2)];
    eta_edge = [t(elem, 3), t(elem, 1)];
    xe_edge = sort(xe_edge, 2);
    eta_edge = sort(eta_edge, 2);
    
    in_flag = 0; out_flag = 0;
    for edge = 1:length(edges(:,1))
        
        for in = 1:length(ein(:, 1))     % is it on the ein boundary?
            if edges(edge, 1) == ein(in, 1) && edges(edge, 2) == ein(in, 2)
                in_flag = 1;
                for l = 1:length(wt1)
                    if edges(edge, 1) == xe_edge(1) && edges(edge, 2) == xe_edge(2)
                        xe = qp1(l); eta = 0;
                    else
                        xe = 0; eta = qp1(l);
                    end
                    [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(xe, eta, num_nodes_per_elem, p, t, elem);
                    dx = p(edges(edge, 1), 1) - p(edges(edge, 2), 1);
                    dy = p(edges(edge, 1), 2) - p(edges(edge, 2), 2);
                    J = sqrt(dx^2 + dy^2);
                    bin = bin + wt1(l) * N * transpose(N) * J;
                    bright = bright + wt1(l) * N * J;
                end
            end
        end
        
        % is it on the eout boundary?
        for in = 1:length(eout(:, 1))
            if edges(edge, 1) == eout(in, 1) && edges(edge, 2) == eout(in, 2)
                out_flag = 1;
                for l = 1:length(wt1)
                    if edges(edge, 1) == xe_edge(1) && edges(edge, 2) == xe_edge(2)
                        xe = qp1(l); eta = 0;
                    else
                        xe = 0; eta = qp1(l);
                    end
                    [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(xe, eta, num_nodes_per_elem, p, t, elem);
                    dx = p(edges(edge, 1), 1) - p(edges(edge, 2), 1);
                    dy = p(edges(edge, 1), 2) - p(edges(edge, 2), 2);
                    J = sqrt(dx^2 + dy^2);
                    bout = bout + wt1(l) * N * transpose(N) * J;
                end
            end
        end
        
        % is it on the Wall boundary? - do nothing, homogeneous Neumann
    end
    
    % place the elemental k matrix into the global K matrix
    for mm = 1:length(perm(:,1))
       i = perm(mm,1);
       j = perm(mm,2);
       K(t(elem, i), t(elem, j)) = K(t(elem, i), t(elem, j)) + k(i,j);
       M(t(elem, i), t(elem, j)) = M(t(elem, i), t(elem, j)) + m(i,j);
       if (in_flag)  
           Bin(t(elem, i), t(elem, j)) = Bin(t(elem, i), t(elem, j)) + bin(i,j);
       end
       if (out_flag)
           Bout(t(elem, i), t(elem, j)) = Bout(t(elem, i), t(elem, j)) + bout(i,j);
       end
    end

    % place the elemental f vector into the global F vector
    for i = 1:num_nodes_per_elem
       if (in_flag)
           F(t(elem, i)) = F((t(elem, i))) + bright(i);
       end
    end
end

bin = F;
end