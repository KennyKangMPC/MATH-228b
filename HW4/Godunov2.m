function [F_left, F_right] = Godunov2(left_point, mid_point, right_point, p_max, u_max)
flux_left = flux(left_point, p_max, u_max);
flux_mid = flux(mid_point, p_max, u_max);
flux_right = flux(right_point, p_max, u_max);

if (mid_point < 0.5 && right_point < 0.5) || (mid_point > 0.5 && right_point > 0.5)
    if mid_point < right_point
        F_right = min(flux_mid, flux_right);
    else
        F_right = max(flux_mid, flux_right);
    end
else
    if mid_point < right_point
        F_right = min(flux_mid, flux_right);
    else
        F_right = flux(p_max / 2, p_max, u_max);
    end
end

if (left_point < 0.5 && mid_point < 0.5) || (left_point > 0.5 && mid_point > 0.5)
    if left_point < mid_point
        F_left = min(flux_left, flux_mid);
    else
        F_left = max(flux_left, flux_mid);
    end
else
    if left_point < mid_point
        F_left = min(flux_left, flux_mid);
    else
        F_left = flux(p_max / 2, p_max, u_max);
    end
end

end

