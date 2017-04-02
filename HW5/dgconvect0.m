%DGCONVECT0  1-D Linear Convection, DG and RK4
n = 10; % number of elements
p = 1; % order of the shape functions
T = 1; % end simulation time
dt = 0.1; % time step size

% UC Berkeley Math 228B, Per-Olof Persson <persson@berkeley.edu>

    h = 1/n;
    s = h * (0:p)' / p;                            % Equidistant nodes
    x = s * ones(1,n) + ones(p+1,1) * (0:h:1-h);   % Entire mesh
    uinit = @(x) exp(-(x - 0.5).^2 / 0.1^2);       % Initial solution
    
    switch p
      case 1
        Mel = h/6 * [2,1; 1,2];                    % Elemental mass-matrix
        Kel = [-1 -1; 1 1]/2;                      % Elemental stiffness-matrix
      case 2
        Mel = h/30 * [4,2,-1; 2,16,2; -1,2,4];     % Elemental mass-matrix
        Kel = [-3,-4,1; 4,0,-4; -1,4,3]/6;         % Elemental stiffness-matrix
      otherwise
        error('Polynomial order not implemented.');
    end
    
    % RHS function in semi-discrete system
    function r = rhs(u)
        r = Kel * u;                               % Integral
        r(end,:) = r(end,:) - u(end,:);            % Right-end flux
        r(1,:) = r(1,:) + u(end, [end, 1:end-1]);  % Left-end flux
        r = Mel \ r;                               % Mass matrix inverse
    end

    u = uinit(x);                                  % Interpolate initial condition
    for it = 1:T/dt                                % Main time-stepping loop
        % Runge-Kutta 4
        k1 = dt * rhs(u);
        k2 = dt * rhs(u + k1/2);
        k3 = dt * rhs(u + k2/2);
        k4 = dt * rhs(u + k3);
        u = u + (k1 + 2*k2 + 2*k3 + k4) / 6;

        % Plotting
        if mod(it, T/dt/1000) == 0                 % 1000 frames
            xx = linspace(0, 1, 1000);             % Fine grid for exact solution
            uexact = uinit(mod(xx - dt*it, 1.0));  % Exact solution
            plot(x, u, 'b', xx, uexact, 'k')
            grid on, axis([0, 1, -0.1, 1.1]), drawnow
        end
    end

    uexact = uinit(mod(x - T, 1.0));               % Exact final solution
    error = max(abs(u(:) - uexact(:)));            % Discrete inf-norm error
