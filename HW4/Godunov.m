function [F_left, F_right] = Godunov(flux_left, flux_mid, flux_right, left_point, mid_point, right_point)

if left_point < mid_point
    F_left = min(flux_left, flux_mid);
else
    F_left = max(flux_left, flux_mid);
end

if mid_point < right_point
    F_right = min(flux_mid, flux_right);
else
    F_right = max(flux_mid, flux_right);
end

end

