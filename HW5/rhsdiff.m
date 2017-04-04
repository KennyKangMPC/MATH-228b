function [r] = rhs(u, Kel, Mel, sigma, k)
    r = Kel * (u - k * sigma);                                  % Integral
    r(end, :) = r(end, :) - (u(end, :) - k * sigma(end, :));    % right flux
    r(1, :) = r(1, :) + (u(end, [end, 1:end-1]) - ...
                k * sigma(end, [end, 1:end-1]));                % left flux
    %r(end,:) = r(end,:) - u(end,:);            % Right-end flux
    %r(1,:) = r(1,:) + u(end, [end, 1:end-1]);  % Left-end flux
    r = Mel \ r;
end

