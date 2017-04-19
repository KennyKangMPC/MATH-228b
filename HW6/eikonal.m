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

F   = 1.0 .* ones(numx, numy);
phi = 0.0 .* ones(numx, numy);
initx = 0.2;
inity = 0.2;
norm = 1.0; % arbitrary initial value
tol = 1e-6; % iteration tolerance
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
    if ifig == 10 % plot every 10 frames
        ifig = 0;
        % plot the contours
        %contour(XX, YY, phi)
        % plot the surface
        hSurf = surf(XX,YY,phi,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        drawnow
    end
    
    % apply the only boundary condition
    phi(initx/dx, inity/dy) = 0.0;
    
    % find the L2-norm of the solution
    norm = sum(sum(abs(phi_prev - phi) .^ 2));
end












