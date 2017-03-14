function [F_left, F_right] = Roe(flux_left, flux_mid, flux_right, left_point, mid_point, right_point)
p_max = 1.0;
u_max = 1.0;

F_left = 0.5 * (flux_left + flux_mid) - 0.5 * u_max * abs(1 - (left_point + mid_point)/p_max) * (mid_point - left_point);
F_right = 0.5 * (flux_mid + flux_right) - 0.5 * u_max * abs(1 - (mid_point + right_point)/p_max) * (right_point - mid_point);

end

