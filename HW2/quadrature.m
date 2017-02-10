function [wt, qp] = quadrature(qp_order)

switch qp_order
    case 1
        wt = [0.5];
        qp = [1/6];
    otherwise
        disp('You entered an unsupported shape function order for the quadrature rule.');
end

