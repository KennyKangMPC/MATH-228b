function [N, dN_dxe, dN_deta, x_xe_eta, y_xe_eta, dx_dxe, dx_deta, dy_dxe, dy_deta, B] = shapefunctions(xe, eta, num_nodes_per_elem, coordinates, LM, elem)

N = zeros(num_nodes_per_elem, 1);
dN = zeros(num_nodes_per_elem, 1);

switch num_nodes_per_elem
    case 3
        N(1) = 1 - xe - eta;
        N(2) = xe;
        N(3) = eta;
        dN_dxe(1) = -1;
        dN_dxe(2) = 1;
        dN_dxe(3) = 0;
        dN_deta(1) = -1;
        dN_deta(2) = 0;
        dN_deta(3) = 1;
        B = [dN_dxe(1), dN_dxe(2), dN_dxe(3); dN_deta(1), dN_deta(2), dN_deta(3)];
    case 6
        N(1) = 1 - 3*xe - 3*eta + 2*xe^2 + 2*eta^2 + 4*xe*eta;
        N(2) = 4*xe - 4*xe^2 - 4*xe*eta;
        N(3) = -xe + 2*xe^2;
        N(4) = 4*xe*eta;
        N(5) = -eta + 2*eta^2;
        N(6) = 4*eta - 4*eta^2 - 4*xe*eta;
        dN_dxe(1) = -3 + 4*xe + 4*eta;
        dN_dxe(2) = 4 - 8*xe - 4*eta;
        dN_dxe(3) = -1 + 4*xe;
        dN_dxe(4) = 4*eta;
        dN_dxe(5) = 0;
        dN_dxe(6) = -4*eta;
        dN_deta(1) = -3 + 4*eta + 4*xe;
        dN_deta(2) = -4*xe;
        dN_deta(3) = 0;
        dN_deta(4) = 4*xe;
        dN_deta(5) = -1 + 4*eta;
        dN_deta(6) = 4 - 8*eta - 4*xe;

%           N(1) = 1 - 3*xe - 3*eta + 2*xe^2 + 2*eta^2 + 4*xe*eta;
%           N(2) = -xe + 2*xe^2;
%           N(3) = -eta + 2*eta^2;
%           N(4) = 4*eta - 4*eta^2 - 4*xe*eta;
%           N(5) = 4*xe - 4*xe^2 - 4*xe*eta;
%           N(6) = 4*xe*eta;
%           dN_dxe(1) = -3 + 4*xe + 4*eta;
%           dN_dxe(2) = -1 + 4*xe;
%           dN_dxe(3) = 0;
%           dN_dxe(4) = -4*eta;
%           dN_dxe(5) = 4 - 8*xe - 4*eta;
%           dN_dxe(6) = 4*eta;
%           dN_deta(1) = -3 + 4*eta + 4*xe;
%           dN_deta(2) = 0;
%           dN_deta(3) = -1 + 4*eta;
%           dN_deta(4) = 4 - 8*eta - 4*xe;
%           dN_deta(5) = -4*xe;
%           dN_deta(6) = 4*xe;

%         N(1) = 1 - 3*xe - 3*eta + 2*xe^2 + 2*eta^2 + 4*xe*eta;
%         N(2) = -eta + 2*eta^2;
%         N(3) = 
%         N(4) = 
%         N(5) = 
%         N(6) = 
        B = [dN_dxe(1), dN_dxe(2), dN_dxe(3), dN_dxe(4), dN_dxe(5), dN_dxe(6); dN_deta(1), dN_deta(2), dN_deta(3), dN_deta(4), dN_deta(5), dN_deta(6)];
    otherwise
        disp('You entered an unsupported number of nodes per element.');
end

% check that all of the shape functions sum to 1 at any location
sum = 0; sum_xe = 0; sum_eta = 0;
for i = 1:num_nodes_per_elem
    sum = sum + N(i);
    sum_xe = sum_xe + dN_dxe(i);
    sum_eta = sum_eta + dN_deta(i);
end

tol = 1e-14;
if abs(sum - 1) > tol || abs(sum_xe - 0) >= tol || abs(sum_eta - 0) >= tol
    disp('Error in shape functions!')
    sprintf('Sum: %f, dxe: %f, deta: %f', sum, sum_xe, sum_eta)
end

% x(xe, eta) and y(xe, eta) transformation to the parametric domain
x_xe_eta = 0.0;
y_xe_eta = 0.0;
dx_dxe = 0.0;
dy_dxe = 0.0;
dx_deta = 0.0;
dy_deta = 0.0;

for i = 1:num_nodes_per_elem
    x_xe_eta = x_xe_eta + coordinates(LM(elem, i), 1) * N(i);
    y_xe_eta = y_xe_eta + coordinates(LM(elem, i), 2) * N(i);
    
    dx_dxe   = dx_dxe   + coordinates(LM(elem, i), 1) * dN_dxe(i);
    dy_dxe   = dy_dxe   + coordinates(LM(elem, i), 2) * dN_dxe(i);
    
    dx_deta = dx_deta   + coordinates(LM(elem, i), 1) * dN_deta(i);
    dy_deta = dy_deta   + coordinates(LM(elem, i), 2) * dN_deta(i);
end