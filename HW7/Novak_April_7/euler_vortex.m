function [r, ru, rv, rE] = euler_vortex(x, y, time, pars)
%EULER_VORTEX  Exact solution to Euler's equations of gas dynamics

% UC Berkeley Math 228B, Per-Olof Persson <persson@berkeley.edu>

gamma = 1.4;
rc    = pars(1); eps   = pars(2); M0    = pars(3);
theta = pars(4); x0    = pars(5); y0    = pars(6);

rinf = 1;
uinf = 1;
Einf = 1/gamma/M0^2/(gamma - 1) + 1/2;
pinf = (gamma - 1) * (Einf - 1/2);
ubar = uinf * cos(theta);
vbar = uinf * sin(theta);
f = (1 - ((x - x0) - ubar*time).^2 - ((y - y0) - vbar*time).^2) ./ rc^2;

u = uinf * (cos(theta) - eps*((y - y0)-vbar*time) / (2*pi*rc) .* exp(f/2));
v = uinf * (sin(theta) + eps*((x - x0)-ubar*time) / (2*pi*rc) .* exp(f/2));
r = rinf * (1 - eps^2 * (gamma - 1) * M0^2/(8*pi^2) .* exp(f)).^(1/(gamma-1));
p = pinf * (1 - eps^2 * (gamma - 1) * M0^2/(8*pi^2) .* exp(f)).^(gamma/(gamma-1));

ru = r.*u;
rv = r.*v;
rE = p/(gamma - 1) + 1/2 * (ru.^2 + rv.^2)./r;
