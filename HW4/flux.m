function [f] = flux(p)
p_max = 1.0;
u_max = 1.0;

% computes the flux for the car example
f = p * u_max * (1 - p/p_max);

end

