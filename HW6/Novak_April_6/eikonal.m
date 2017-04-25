% solve the Eikonal equation on unit square
clear all
X = 1;
Y = 1;

dx = 0.01;  dy = 0.01;      % mesh spacing = 1/100
dt = 1e-3;                  % time step
x  = 0:dx:X; y = 0:dy:Y;    % Cartesian grid
[XX, YY] = meshgrid(x, y);
numx = length(x); 
numy = length(y); 

opt   = 4;                  % flag to select which case to run
initx = 0.2;                % x-coordinate of starting point
inity = 0.2;                % y-coordinate of starting point
norm = 1.0; ifig = 0;       % arbitrary initial values
tol  = 1e-10;               % iteration tolerance

phi = 0.0 .* ones(numx, numy);
F   = ones(numx, numy);
nx  = zeros(numx, numy);
ny  = zeros(numx, numy);

switch opt
    case 1
        arrivx = 0.8;
        arrivy = 0.8;
    case 2
        for i = 1:numy
            if y(i) < 0.5
                F(i, :) = 0.5;
            end
        end
        arrivx = 0.8;
        arrivy = 0.45;
    case 3
        for i = 1:numx
            for j = 1:numy
                F(i, j) = 1 - 0.9 * cos(4 .* pi .* x(j)) ...
                    .* exp(-10 .* ((x(j) - 0.5) .^ 2 + (y(i) - 0.5).^2));
            end
        end
        arrivx = 0.8;
        arrivy = 0.8;
    case 4
        initx = 0.05;
        inity = 0.05;
        arrivx = 0.8;
        arrivy = 0.8;
        low = [0.1, 0.3];
        for i = 1:numx
            for j = 1:numy
                if (low(1) < x(i)) && (x(i) < low(2))
                    if (low(1) < y(j)) && (y(j) < low(2))
                        F(i, j) = 0.01;
                    end
                end
                if (low(1) + 0.4 < x(i)) && (x(i) < low(2) + 0.4)
                    if (low(1) + 0.4 < y(j)) && (y(j) < low(2) + 0.4)
                        F(i, j) = 0.01;
                    end
                end
            end
        end
    otherwise
        disp('Problem selection undefined.')
end
       
while norm > tol
    phi_prev = phi;
    
    for i = 1:length(x) % loop over x-coordinates
        for j = 1:length(y) % loop over y-coordinates
            plus = grad_plus(phi, i, j, dx, dy, numx, numy);
            phi(i, j) = phi(i, j) - dt * (max(F(i, j), 0) * plus) + dt;
        end
    end
    
    % apply the only boundary condition
    phi(initx/dx, inity/dy) = 0.0;
    
    % find the infinity norm of the difference
    norm = max(max(abs(phi_prev - phi)));
end

contour(XX, YY, phi)
hold on

% compute gradients, then find the normal vectors
[fx, fy] = gradient(phi, dx, dy);
nx       = fx ./ sqrt(fx .^ 2 + fy .^ 2);
ny       = fy ./ sqrt(fx .^ 2 + fy .^ 2);

dn   = 0.01; % time step for tracing out the shortest path
xc   = [x(arrivx/dx)];
yc   = [y(arrivy/dy)];
dist = sqrt((xc - initx).^2 + (yc - inity).^2);

% find the optimal path by moving in direction opposite the gradient
i = 1;
while ((dist > 2*dx) && (i < 2000))
    nx_loc = interp2(XX, YY, nx, xc(end), yc(end));
    ny_loc = interp2(XX, YY, ny, xc(end), yc(end));

    i = i + 1;
    if i > 2000
        break;
    end
    
    xc = [xc, xc(end) - nx_loc * dn];
    yc = [yc, yc(end) - ny_loc * dn];
    
    dist = sqrt((xc(end) - initx).^2 + (yc(end) - inity).^2);
end

hold on
scatter(xc, yc, 'ko')

saveas(gcf, sprintf('case%i.png',opt))