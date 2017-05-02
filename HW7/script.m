clear all
N = [32, 64, 128];
inf_norm = zeros(length(N), 4);

for i = 1:length(N)
    n = N(i);
    h = 10 / N(i);
    x = h:h:10; y = h:h:10;
    [X, Y] = meshgrid(x, y);

    final = sqrt(20^2 + 10^2);
    dt = final/(ceil(final / (0.2 * h)));
    dt = dt / 4;
    num_dt = final/dt;

    % obtain the initial condition from euler_vortex
    pars = [0.5, 1, 0.5, atan2(1,2), 5, 5];
    [r, ru, rv, rE] = euler_vortex(X, Y, 0, pars);
    r_ex = r; ru_ex = ru; rv_ex = rv; rE_ex = rE;

    % RK-4 timestepping
    for it = 1:num_dt
        [r, ru, rv, rE] = euler_rk4step(r, ru, rv, rE, h, dt);
        if (mod(it, 100) == 0)
            surf(X, Y, r)
        end
    end

    inf_norm(i, :) = [max(max(abs(r_ex - r))), max(max(abs(ru_ex - ru))), ...
        max(max(abs(rv_ex - rv))), max(max(abs(rE_ex - rE)))];
end

loglog(10./N, inf_norm(:, 1))
hold on
loglog(10./N, inf_norm(:, 2))
hold on
loglog(10./N, inf_norm(:, 3))
hold on
loglog(10./N, inf_norm(:, 4))
hold on
legend('density', 'x-momentum', 'y-momentum', 'total energy')

