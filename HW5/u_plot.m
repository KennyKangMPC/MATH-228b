function [u_plot] = u_plot(n, p, u, x_plot, x)
% determines the solution in the physical domain for plotting

for el = 1:n
    for i = 1:(p + 1)
        A(i, 1:(p + 1)) = x(i, el) .^ ((1:(p + 1)) - 1);
    end

    coeff_plot(:, el) = A \ u(:, el);

    for i = 1:length(x_plot(:, 1))
        u_plot(i, el) = sum(coeff_plot(:, el)' .* (x_plot(i, el) .^ ((1:(p + 1)) - 1)));
    end
end

end

