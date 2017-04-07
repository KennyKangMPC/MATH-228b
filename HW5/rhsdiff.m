function [r] = rhs(u, Kel, Mel, k)
    % solve for sigma
    r = - Kel * u;
    r(end, :) = r(end, :) + u(1, [2:end, 1]);   % right-end flux
    r(1, :) = r(1, :) - u(1, :);                % left-end flux
    sigma = Mel \ r;

    % solve for u
    r = Kel * (u - k * sigma);                                  % Integral
    r(end, :) = r(end, :) - (u(end, :) - k * sigma(end, :));    % right flux
    r(1, :) = r(1, :) + (u(end, [end, 1:end-1]) - ...
                k * sigma(end, [end, 1:end-1]));                % left flux
    r = Mel \ r;
end

