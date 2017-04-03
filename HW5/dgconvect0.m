%DGCONVECT0  1-D Linear Convection, DG and RK4
clear all
n = 10; % number of elements
p = 2; % order of the shape functions
T = 0.25; % end simulation time
dt = 1e-4; % time step size

% UC Berkeley Math 228B, Per-Olof Persson <persson@berkeley.edu>

    h = 1/n;
    s_master = - cos(pi * (0:p) / p);
    s = ((s_master + 1).* h / 2)';                        % Chebyshev nodes
    x = s * ones(1,n) + ones(p+1,1) * (0:h:1-h);   % Entire mesh
    uinit = @(x) exp(-(x - 0.5).^2 / 0.1^2);       % Initial solution
    
    % find basis function coefficients (c_i^j) for each element
    A = []; 
    I = eye(length(s));
    
    % coefficients are found in the master domain
    for i = 1:length(s)
        A(i, :) = legendre_poly(s_master(i), p);
    end
    
    c = I \ A;
    
    % quadrature points and weights
    [qp, wt] = gauss_quad(2*p);
    
    M = zeros(p + 1, p + 1);
    K = zeros(p + 1, p + 1);
    % loop over the quadrature points
    for q = 1:length(qp)
        [yi, dyi] = legendre_poly(qp(q), p);
        [yj, dyj] = legendre_poly(qp(q), p);
        for i = 1:(p + 1)
            for j = 1:(p + 1)
                M(i, j) = M(i, j) + wt(q) * dot(c(i, :), yi) * dot(c(j, :), yj) * h / 2;
                K(i, j) = K(i, j) + wt(q) * dot(c(j, :), yj) * dot(c(i, :), dyi);
            end
        end
    end
    
    Mel = M;
    Kel = K;
    
    u = uinit(x);                                  % Interpolate initial condition
    
    % sum of phi_i(x_k) from i = 0 to P
    [y_right] = legendre_poly(1, p)';
    [y_left] = legendre_poly(-1, p)';
    
    %y_right = y_right * ones(1, n);
    %y_left = y_left * ones(1, n);
    
    %y_right = y_right * u(end, :);
    %y_left = y_left * u(end, [end, 1:(end-1)]);


    
    
    % temporary
    yr = y_right;
    yl = y_left;
    
%     for it = 1:T/dt                                % Main time-stepping loop
%         % Runge-Kutta 4
%         k1 = dt * rhs(u, Kel, Mel, yr, yl);
%         k2 = dt * rhs(u + k1/2, Kel, Mel, yr, yl);
%         k3 = dt * rhs(u + k2/2, Kel, Mel, yr, yl);
%         k4 = dt * rhs(u + k3, Kel, Mel, yr, yl);
%         u = u + (k1 + 2*k2 + 2*k3 + k4) / 6;
% 
%         % Plotting
%         if mod(it, T/dt/1000) == 0                 % 1000 frames
%             xx = linspace(0, 1, 1000);             % Fine grid for exact solution
%             uexact = uinit(mod(xx - dt*it, 1.0));  % Exact solution
%             plot(x, u, 'b', xx, uexact, 'k')
%             grid on
%             %axis([0, 1, -0.1, 1.1])
%             drawnow
%         end
%     end
% 
%     uexact = uinit(mod(x - T, 1.0));               % Exact final solution
%     error = max(abs(u(:) - uexact(:)));            % Discrete inf-norm error
%     
%     % compute the error as the (integrated) L-2 norm
%     
