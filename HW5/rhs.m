function [r] = rhs(u, Kel, Mel)
    r = Kel * u;                                        % Integral
    r(end,:) = r(end,:) - u(end,:);            % Right-end flux
    r(1,:) = r(1,:) + u(end, [end, 1:end-1]);  % Left-end flux
    r = Mel \ r;
end

