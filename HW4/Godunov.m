function [F_left, F_right] = Godunov(left_point, mid_point, right_point, p_max, u_max)
maxpt = p_max/2;
maxflux = flux(maxpt, p_max, u_max);

flux_left = flux(left_point, p_max, u_max);
flux_mid = flux(mid_point, p_max, u_max);
flux_right = flux(right_point, p_max, u_max);

if left_point < mid_point
    if ((maxpt > left_point) && (maxpt < mid_point))
        F_left = min([flux_left, flux_mid, maxflux]);
    else
        F_left = min([flux_left, flux_mid]);
    end
else
    if ((maxpt < left_point) && (maxpt > mid_point))
        F_left = max([flux_left, flux_mid, maxflux]);
    else
        F_left = max([flux_left, flux_mid]);
    end
end

if mid_point < right_point
    if ((maxpt > mid_point) && (maxpt < right_point))
        F_right = min([flux_left, flux_mid, maxflux]);
    else
        F_right = min([flux_mid, flux_right]);
    end
else
    if ((maxpt < mid_point) && (maxpt > right_point))
        F_right = max([flux_left, flux_mid, maxflux]);
    else
        F_right = max([flux_mid, flux_right]);
    end
end

end

