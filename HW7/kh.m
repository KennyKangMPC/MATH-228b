clear all
N = 64;
h = 1 / N;
dt = 0.2 * h;

x = h:h:1; y = h:h:1;
[X, Y] = meshgrid(x, y);

num_dt = 2.0/dt; 
P = 3; gamma = 7/5;

% set the initial condition
r = zeros(length(x), length(y));
for i = 1:length(y)
    for j = 1:length(x)
        if (abs(y(i) - 0.5) < (0.15 + sin(2 * pi * x(j)) / 200))
            r(i, j) = 2;
        else
            r(i, j) = 1;
        end
    end
end

ru = r .* (r - 1);
rv = 0 .* r;
rE = P/(gamma - 1) + (ru.^2 + rv.^2)./(2.*r);

% % RK-4 timestepping
% for it = 1:num_dt
%     [r, ru, rv, rE] = euler_rk4step(r, ru, rv, rE, h, dt);
%     if (mod(it, 100) == 0)
%         surf(X, Y, r)
%     end
% end