% solve the Eikonal equation on unit square
clear all
X = 1;
Y = 1;

dx = 0.01;  dy = 0.01;      % mesh spacing = 1/100
dt = 1e-4;                  % time step
x  = 0:dx:X; y = 0:dy:Y;    % Cartesian grid
[XX, YY] = meshgrid(x, y);
numx = length(x); 
numy = length(y); 

opt   = 3;                  % flag to select which case to run
plot  = 1;                  % flag to select real-time plot
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
                F(i, j) = 1 - 0.9 * cos(4 .* pi .* x(i)) ...
                    .* exp(-10 .* ((x(i) - 0.5) .^ 2 + (y(j) - 0.5).^2));
            end
        end
        arrivx = 0.8;
        arrivy = 0.8;
    case 4
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
            end
        end
        
        % round out the corners of the box
        F(low(1)/dx, low(1)/dy) = 1.0;
        F(low(1)/dx, low(2)/dy) = 1.0;
        F(low(2)/dx, low(1)/dy) = 1.0;
        F(low(2)/dx, low(2)/dy) = 1.0;
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
    
    ifig = ifig + 1;
    if (mod(ifig, 20) == 0) && (plot) % plot every 20 frames
        % plot the contours
        contour(XX, YY, phi)
        % plot the surface
        %hSurf = surf(XX,YY,phi,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        drawnow
    end
    
    % apply the only boundary condition
    phi(initx/dx, inity/dy) = 0.0;
    
    % find the L2-norm of the solution
    norm = sum(sum(abs(phi_prev - phi) .^ 2));
end

if (~plot)
    contour(XX, YY, phi)
end

% compute gradients, then find the normal vectors
[fx, fy] = gradient(phi, dx, dy);
nx       = fx ./ sqrt(fx .^ 2 + fy .^ 2);
ny       = fy ./ sqrt(fx .^ 2 + fy .^ 2);

hold on
% only plot a subsection of the normals so it's easier to see
quiver(XX(1:2:end, 1:2:end), YY(1:2:end, 1:2:end), ...
    nx(1:2:end, 1:2:end), ny(1:2:end, 1:2:end))

dn   = 0.01; % time step for tracing out the shortest path
xc   = [x(arrivx/dx)];
yc   = [y(arrivy/dy)];
dist = sqrt((xc - arrivx).^2 + (yc - arrivy).^2);
i    = 0;

% find the optimal path
while ((dist > dx) && (i < 1000))
    nx_loc = interp2(XX, YY, nx, xc(end), yc(end));
    ny_loc = interp2(XX, YY, ny, xc(end), yc(end));

    % move in the direction opposite the gradient
    xc = [xc, xc(end) - nx_loc * dn];
    yc = [yc, yc(end) - ny_loc * dn];
    
    i = i + 1;
    %if (mod(i, 1) == 0) % only plot every x points
    %    plot(xc, yc, 'ko')
    %    hold on
    %end
    
    dist = sqrt((xc(end) - initx).^2 + (yc(end) - inity).^2);
end

scatter(xc, yc, 'ko')

saveas(gcf, sprintf('case%i.png',opt))