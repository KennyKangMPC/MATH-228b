% Q1, HW 4, using Godunov's method
clear all
p_max = 1.0;
u_max = 1.0;
p_left = 0.8;
p_right = 0.0;
dx = 4/400;
dt = 0.8 * dx / u_max;
X = 5;
T = 3;

mesh = 0:dx:X;
time = 0:dt:T;

p = cell(length(time), 1);
for i = 1:length(time)
    p{i} = p_right .* ones(size(mesh));
end

for t = 1:length(time)
    F_left = zeros(size(mesh));
    F_right = zeros(size(mesh));

    % sweep over all of the nodes for a single time step, except the last node
    % and compute the left and right fluxes
    for i = 1:(length(mesh) - 1)
        if i == 1
            left_point = p_left;
        else
            left_point = p{t}(i - 1);
        end
        
        flux_left = flux(left_point);
        flux_mid = flux(p{t}(i));
        flux_right = flux(p{t}(i + 1));
        
        [F_left(i), F_right(i)] = Godunov(flux_left, flux_mid, flux_right, left_point, p{t}(i), p{t}(i+1));

    p{t + 1}(i) = p{t}(i) - (dt / dx) * (F_right(i) - F_left(i));
    end
end

for t = 1:length(time)
    plot(mesh, p{t});
    drawnow
    pause(0.00001);
end
