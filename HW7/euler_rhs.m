function [fr, fru, frv, frE] = euler_rhs(r, ru, rv, rE, h)

% determine the Euler fluxes
[Frx, Fry, Frux, Fruy, Frvx, Frvy, FrEx, FrEy] = euler_fluxes(r, ru, rv, rE);

% determine the spectral divergence
fr = spectral_divergence(Frx, Fry, h);
fru = spectral_divergence(Frux, Fruy, h);
frv = spectral_divergence(Frvx, Frvy, h);
frE = spectral_divergence(FrEx, FrEy, h);