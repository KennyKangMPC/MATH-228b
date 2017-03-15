function [F_left, F_right] = Roe(left_point, mid_point, right_point, p_max, u_max)
flux_left = flux(left_point, p_max, u_max);
flux_mid = flux(mid_point, p_max, u_max);
flux_right = flux(right_point, p_max, u_max);

F_left = 0.5 * (flux_left + flux_mid) - ...
    0.5 * u_max * abs(1 - (left_point + mid_point)/p_max) * (mid_point - left_point);
F_right = 0.5 * (flux_mid + flux_right) - ...
    0.5 * u_max * abs(1 - (mid_point + right_point)/p_max) * (right_point - mid_point);
end

