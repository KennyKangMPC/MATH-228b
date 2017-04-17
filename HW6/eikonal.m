% solve the Eikonal equation on unit square
X = 1;
Y = 1;

dx = 0.1;           % mesh spacing
dy = 0.1;
dt = 0.1;           % time step
x = 0:dx:X;         % Cartesian grid
y = 0:dy:Y;         % Cartesian grid
numx = length(x);   
numy = length(y);   

F   = 1.0 .* ones(numx, numy);
phi = zeros(numx, numy);

for i = 1:length(x) % loop over x-coordinates
    for j = 1:length(y) % loop over y-coordinates
        plus = grad_plus(phi, i, j, dx, dy, numx, numy);
        phi(i, j) = phi(i, j) - dt * (max(F(i, j), 0) * plus + min(F(i, j), 0)) + dt;
    end
end











