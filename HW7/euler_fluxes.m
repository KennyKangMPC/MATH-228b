function [Frx, Fry, Frux, Fruy, Frvx, Frvy, FrEx, FrEy] = euler_fluxes(r, ru, rv, rE)

u = ru ./ r; v = rv ./ r; E = rE ./ r;

% compute pressure from ideal gas law
y = 7/5;
P = (y - 1) .* r .* (E - (u.^2 + v.^2) / 2);

Frx = ru;               Fry = rv;
Frux = r .* u.^2 + P;   Fruy = ru .* v;
Frvx = ru .* v;         Frvy = r .* v.^2 + P;
FrEx = u .* (rE + P);   FrEy = v .* (rE + P);
