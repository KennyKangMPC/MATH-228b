function [r, ru, rv, rE] = euler_rk4step(r, ru, rv, rE, h, dt)

% determine k1, k2, k3, k4
[k1_fr, k1_fru, k1_frv, k1_frE] = euler_rhs(r, ru, rv, rE, h);
k1_fr  = dt .* k1_fr;
k1_fru = dt .* k1_fru;
k1_frv = dt .* k1_frv;
k1_frE = dt .* k1_frE;

[k2_fr, k2_fru, k2_frv, k2_frE] = euler_rhs(r + k1_fr./2, ru + k1_fru./2, rv + k1_frv./2, rE + k1_frE./2, h);
k2_fr  = dt .* k2_fr;
k2_fru = dt .* k2_fru;
k2_frv = dt .* k2_frv;
k2_frE = dt .* k2_frE;

[k3_fr, k3_fru, k3_frv, k3_frE] = euler_rhs(r + k2_fr./2, ru + k2_fru./2, rv + k2_frv./2, rE + k2_frE./2, h);
k3_fr  = dt .* k3_fr;
k3_fru = dt .* k3_fru;
k3_frv = dt .* k3_frv;
k3_frE = dt .* k3_frE;

[k4_fr, k4_fru, k4_frv, k4_frE] = euler_rhs(r + k3_fr, ru + k3_fru, rv + k3_frv, rE + k3_frE, h);
k3_fr  = dt .* k3_fr;
k3_fru = dt .* k3_fru;
k3_frv = dt .* k3_frv;
k3_frE = dt .* k3_frE;

r  =  r + (k1_fr + 2 .* k2_fr + 2 .* k3_fr + k4_fr) ./ 6;
ru = ru + (k1_fru + 2 .* k2_fru + 2 .* k3_fru + k4_fru) ./ 6;
rv = rv + (k1_frv + 2 .* k2_frv + 2 .* k3_frv + k4_frv) ./ 6;
rE = rE + (k1_frE + 2 .* k2_frE + 2 .* k3_frE + k4_frE) ./ 6;

% filter each solution component (x and y for the four parts)
r = spectral_filter(r);
ru = spectral_filter(ru);
rv = spectral_filter(rv);
rE = spectral_filter(rE);

end

