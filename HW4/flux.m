function [f] = flux(p, p_max, u_max)
f = p .* u_max .* (1 - p./p_max);
end

