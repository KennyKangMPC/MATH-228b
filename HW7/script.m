n = 64;
h = 10/n;

x = h:h:10; % square domain, 0 < x, y < 10
y = h:h:10;
[X, Y] = meshgrid(x, y);

final = sqrt(20^2 + 10^2);
dt = final/(ceil(final / (0.2 * h)));
dt = dt / 2;
num_dt = final/dt;

% obtain the initial condition from euler_vortex
pars = [0.5, 1, 0.5, atan2(1,2), 5, 5];
[r, ru, rv, rE] = euler_vortex(X, Y, 0, pars);

% RK-4 timestepping
for it = 1:num_dt
    [r, ru, rv, rE] = euler_rk4step(r, ru, rv, rE, h, dt);
    [r_ex, ru_ex, rv_ex, rE_ex] = euler_vortex(X, Y, it * dt, pars);
    surf(X, Y, r_ex)
    drawnow
end

% determine the exact final solution
[r_ex, ru_ex, rv_ex, rE_ex] = euler_vortex(X, Y, 0, pars);

inf_norm = max([max(max(r_ex - r)), max(max(ru_ex - ru)), max(max(rv_ex - rv)), max(max(rE_ex - rE))]);


