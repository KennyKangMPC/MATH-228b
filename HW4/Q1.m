% Q1, HW 4
clear all
p_max = 1.0;
u_max = 1.0;
p_left = 0.8;
p_right = 0.0;
dx = 4/400;
dt = 0.8 * dx / u_max;

% problem domain - node locations
mesh = 0:dx:5;

p = zeros(size(mesh));
F_left = zeros(size(mesh));
F_right = zeros(size(mesh));
i = 1;

% compute F_{i - 1/2} for the very first grid point - requires ghost cells
flux_left = flux(p_left);
flux_mid = flux(p(i));
flux_right = flux(p(i + 1));

if p_left < p(i)
    F_left(i) = min(flux_left, flux_mid);
else
    F_left(i) = max(flux_left, flux_mid);
end

if p(i) < p(i+1)
    F_right(i) = min(flux_mid, flux_right);
else
    F_right(i) = max(flux_mid, flux_right);
end

% sweep over all of the nodes for a single time step, except the last node
for i = 2:(length(mesh) - 1)
    flux_left = flux(p(i - 1));
    flux_mid = flux(p(i));
    flux_right = flux(p(i + 1));

    if p(i - 1) < p(i)
        F_left(i) = min(flux_left, flux_mid);
    else
        F_left(i) = max(flux_left, flux_mid);
    end

    if p(i) < p(i + 1)
        F_right(i) = min(flux_mid, flux_right);
    else
        F_right(i) = max(flux_mid, flux_right);
    end
    
end

