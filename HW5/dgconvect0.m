%DGCONVECT0  1-D Linear Convection, DG and RK4
clear all
n = 80; % number of elements
p = 1; % order of the shape functions
T = 0.1; % end simulation time
dt = 1e-4; % time step size
fine_el = 100; % number of points to resolve per discontinuous element
exact_el = fine_el * n - (n - 1);

h = 1/n;
s_master = - cos(pi * (0:p) / p);
s = ((s_master + 1).* h / 2)';                 % Chebyshev nodes
x = s * ones(1,n) + ones(p+1,1) * (0:h:1-h);   % Entire mesh
uinit = @(x) exp(-(x - 0.5).^2 / 0.1^2);       % Initial solution
xx = linspace(0, 1, exact_el);                     % Fine grid for exact soln

% find basis function coefficients (c_i^j) for each element
A = []; 
I = eye(length(s));

% coefficients are found in the master domain
for i = 1:length(s)
    [y, dy] = legendre_poly(s_master(i), p);
    A(i, :) = y;
end

c = A \ I;

% quadrature points and weights
[qp, wt] = gauss_quad(2*p);

Mel = zeros(p + 1, p + 1);
Kel = zeros(p + 1, p + 1);

for q = 1:length(qp)
    [yi, dyi] = legendre_poly(qp(q), p);
    [yj, dyj] = legendre_poly(qp(q), p);
    for i = 1:(p + 1)
        for j = 1:(p + 1)
            Mel(i, j) = Mel(i, j) + wt(q) * dot(c(:, i), yi) * dot(c(:, j), yj) * h / 2;
            Kel(i, j) = Kel(i, j) + wt(q) * dot(c(:, j), yj) * dot(c(:, i), dyi);
        end
    end
end

u = uinit(x);                                  % Interpolate initial condition

% create matrix of evaluation points for plotting
x_plot = [];
for el = 1:n
    x_plot(:, el) = linspace(x(1, el), x(end, el), 3*p);
    endg

horiz = [[0, x(end, :)]; [0, x(end, :)]];
vert = [0; 1] * ones(1, length([0, x(end, :)]));

for it = 1:T/dt                                % Main time-stepping loop
    % Runge-Kutta 4
    k1 = dt * rhs(u, Kel, Mel);
    k2 = dt * rhs(u + k1/2, Kel, Mel);
    k3 = dt * rhs(u + k2/2, Kel, Mel);
    k4 = dt * rhs(u + k3, Kel, Mel);
    u = u + (k1 + 2*k2 + 2*k3 + k4) / 6;

    % Plotting
    if mod(it, T/dt/500) == 0                 % 1000 frames

        uexact = uinit(mod(xx - dt*it, 1.0));  % Exact solution
        u_plot = []; A = []; coeff_plot = [];

        for el = 1:n
            for i = 1:(p + 1)
                A(i, 1:(p + 1)) = x(i, el) .^ ((1:(p + 1)) - 1);
            end

            coeff_plot(:, el) = A \ u(:, el);

            for i = 1:length(x_plot(:, 1))
                u_plot(i, el) = sum(coeff_plot(:, el)' .* (x_plot(i, el) .^ ((1:(p + 1)) - 1)));
            end
        end

        %plot(x_plot, u_plot, 'b*-', x, u, 'r', xx, uexact, 'k', horiz, vert, '--')
        %grid on
        %axis([0, 1, -0.1, 1.1])
        %drawnow
    end
end

uexact = uinit(mod(x - T, 1.0));               % Exact final solution
error = max(abs(u(:) - uexact(:)));            % Discrete inf-norm error
   



% determines the solution in the physical domain for plotting
x_norm = [];
for el = 1:n
    x_norm(:, el) = linspace(x(1, el), x(end, el), fine_el);
end

A = []; coeff_plot = []; u_norm = [];
for el = 1:n
    for i = 1:(p + 1)
        A(i, 1:(p + 1)) = x(i, el) .^ ((1:(p + 1)) - 1);
    end

    coeff_plot(:, el) = A \ u(:, el);

    for i = 1:length(x_norm(:, 1))
        u_norm(i, el) = sum(coeff_plot(:, el)' .* (x_norm(i, el) .^ ((1:(p + 1)) - 1)));
    end
end

% u_norm : DG solution
% u_exact_norm : exact solution
u_exact_n = uinit(mod(xx - T, 1.0));
u_exact_norm = [];

i = 1;
j = fine_el;
for el = 1:n
    u_exact_norm(:, el) = u_exact_n(i:j)';
    i = i + fine_el - 1;
    j = j + fine_el - 1;
end

% L-2 norm for each element is determined using trapz() at the last time
% step only (follows the original code)
L2_elem = [];
for el = 1:n
    L2_elem(el) = trapz(x_norm(:, el), (u_exact_norm(:, el)' - u_norm(:, el)') .^ 2);
end

% sum over all the elements
L2_norm = sqrt(sum(L2_elem))

