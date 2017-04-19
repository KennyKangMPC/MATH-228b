% solve the Eikonal equation on unit square
clear all
X = 1;
Y = 1;

dx = 0.01;           % mesh spacing
dy = 0.01;
dt = 1e-3;           % time step
x = 0:dx:X;         % Cartesian grid
y = 0:dy:Y;         % Cartesian grid
[XX, YY] = meshgrid(x, y);
numx = length(x);   
numy = length(y); 

opt = 3;
F = ones(numx, numy);
initx = 0.2;
inity = 0.2;

switch opt
    case 1
    case 2
        for i = 1:numy
            if y(i) < 0.5
                F(i, :) = 0.5;
            end
        end
    case 3
        for i = 1:numx
            for j = 1:numy
                F(i, j) = 1 - 0.9 * cos(4 .* pi .* x(i)) ...
                    .* exp(-10 .* ((x(i) - 0.5) .^ 2 + (y(j) - 0.5).^2));
            end
        end
    otherwise
        disp('Problem selection undefined.')
end
        
phi = 0.0 .* ones(numx, numy);
norm = 1.0; % arbitrary initial value
tol = 1e-7; % iteration tolerance
ifig = 0;

while norm > tol
    phi_prev = phi;
    
    for i = 1:length(x) % loop over x-coordinates
        for j = 1:length(y) % loop over y-coordinates
            plus = grad_plus(phi, i, j, dx, dy, numx, numy);
            phi(i, j) = phi(i, j) - dt * (max(F(i, j), 0) * plus) + dt;
        end
    end
    
    ifig = ifig + 1;
    if ifig == 20 % plot every 10 frames
        ifig = 0;
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












