function [F_left, F_right] = Roe(left_point, mid_point, right_point, p_max, u_max)

flux_left = flux(left_point);
flux_mid = flux(mid_point);
flux_right = flux(right_point);

F_left = 0.5 * (flux_left + flux_mid) - 0.5 * u_max * abs(1 - (left_point + mid_point)/p_max) * (mid_point - left_point);
F_right = 0.5 * (flux_mid + flux_right) - 0.5 * u_max * abs(1 - (mid_point + right_point)/p_max) * (right_point - mid_point);

%a_minus = u_max * abs(1 - (left_point + mid_point)/(2 * p_max));
%a_plus = u_max * abs(1 - (right_point + mid_point)/(2 * p_max));

%a_minus = abs(u_max * (flux_mid - flux_left) / (mid_point - left_point));
%a_plus = abs(u_max * (flux_right - flux_mid) / (right_point - mid_point));

%F_left = 0.5 * (flux_left + flux_mid) - 0.5 * abs(a_minus) * (mid_point - left_point);
%F_right = 0.5 * (flux_mid + flux_right) - 0.5 * abs(a_plus) * (right_point - mid_point);

end

