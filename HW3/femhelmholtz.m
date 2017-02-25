
wave = 6;              % wave number
imag = sqrt(-1);

% --- simple test case
%  p = [0,0; 1,0; 2,0; 0,1; 1,1; 2,1; 0,2; 1,2; 2,2];
%  t = [1,2,4; 2,5,4; 2,3,5; 3,6,5; 4,5,7; 5,8,7; 5,6,8; 6,9,8];
%  In = [1, 4; 4, 7];
%  Out = [3, 6; 6, 9];
% ----

if (1)
    % % generate the mesh
    pv = [0,0; 5,0; 5,1; 0,1; 0,0];
    [p, t, e] = pmesh(pv, 0.15, 1);
    [In, Out, Wall] = waveguide_edges(p, t);
end


LM = t;                         % location matrix                
num_elem = length(LM(:,1));     % number of elements
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

% assemble the elemental k and elemental f
K = zeros(num_nodes); M = zeros(num_nodes);
Bin = zeros(num_nodes); Bout = zeros(num_nodes); F = zeros(num_nodes, 1);

% change the LM for weird elements
[LM] = changeLM(num_elem, LM, In, Out);

for elem = 1:num_elem
    k = 0; m = 0; bout = 0; bin = 0; bright = 0;

    for ll = 1:length(wt)
         [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(qp(ll, 1), qp(ll, 2), num_nodes_per_elem, p, LM, elem);
         F_mat = transpose([dx_dxe, dx_deta; dy_dxe, dy_deta]);
         J = det(F_mat);
         D = inv(F_mat) * B;

         k = k + wt(ll) * transpose(D) * (D) * J;
         m = m + wt(ll) * N * transpose(N) * J; 
    end
    
    % find the edges of the current element
    edges = [LM(elem,[1,2]); LM(elem,[2,3]); LM(elem,[3,1])];
    edges = sort(edges, 2);
    
    xe_edge = [LM(elem, 1), LM(elem, 2)];
    eta_edge = [LM(elem, 3), LM(elem, 1)];
    xe_edge = sort(xe_edge, 2);
    eta_edge = sort(eta_edge, 2);
    
    in_flag = 0; out_flag = 0;
    for edge = 1:length(edges(:,1))
        
        for in = 1:length(In(:, 1))     % is it on the In boundary?
            if edges(edge, 1) == In(in, 1) && edges(edge, 2) == In(in, 2)
                in_flag = 1;
                for l = 1:length(wt1)
                    if edges(edge, 1) == xe_edge(1) && edges(edge, 2) == xe_edge(2)
                        xe = qp1(l); eta = 0;
                    else
                        xe = 0; eta = qp1(l);
                    end
                    [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(xe, eta, num_nodes_per_elem, p, LM, elem);
                    dx = p(edges(edge, 1), 1) - p(edges(edge, 2), 1);
                    dy = p(edges(edge, 1), 2) - p(edges(edge, 2), 2);
                    J = sqrt(dx^2 + dy^2);
                    bin = bin + wt1(l) * N * transpose(N) * J;
                    bright = bright + wt1(l) * N * J;
                end
            end
        end
        
        % is it on the Out boundary?
        for in = 1:length(Out(:, 1))
            if edges(edge, 1) == Out(in, 1) && edges(edge, 2) == Out(in, 2)
                out_flag = 1;
                for l = 1:length(wt1)
                    if edges(edge, 1) == xe_edge(1) && edges(edge, 2) == xe_edge(2)
                        xe = qp1(l); eta = 0;
                    else
                        xe = 0; eta = qp1(l);
                    end
                    [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(xe, eta, num_nodes_per_elem, p, LM, elem);
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

Kk = K - (wave.^2) .* M + imag .* wave .* (Bin + Bout);
Ff = F .* 2 .* imag .* wave;

a = Kk\Ff;

% plot the exact solution
real_exact = zeros(length(a), 1);
imag_exact = zeros(length(a), 1);

for i = 1:length(real_exact(:, 1))
    real_exact(i) = cos(wave .* p(i, 1));
    imag_exact(i) = imag .* sin(wave .* p(i, 1));
end

real_fem = real(a);

tplot(p, LM, real_fem - real_exact)