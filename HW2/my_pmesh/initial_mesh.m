function [pv] = initial_mesh(pv, hmax, pv_orig)
% place points on the domain boundaries according to hmax

k = 1;
l = 1;
new_pts = [];

for i = 1:(size(pv(:,1)) - 1)
    p = 1;
    x0 = pv(i, 1);
    x_next = pv(i + 1, 1);
    y0 = pv(i, 2);
    y_next = pv(i + 1, 2);
    
    if x_next == x0         % vertical line
        dist = abs(y_next - y0);
        num_divs = ceil(dist/hmax);
        spacing = dist/num_divs;
        
        for j = 1:(num_divs - 1)
            if y_next > y0
                y = y0 + spacing;
                y0 = y;
            else
                y = y0 - spacing;
                y0 = y;
            end
            new_pts(k, :) = [x0, y];
            k = k + 1;
        end
    else        % not a vertical line
        m(l) = (y_next - y0) / (x_next - x0);
        
        % determine number of points to place to get close to hmax
        dist = sqrt((y_next - y0)^2 + (x_next - x0)^2);
        num_divs = ceil(dist/hmax);
        spacing = dist/num_divs;
        
        % compute the new points
        for j = 1:(num_divs - 1)
            if x_next > x0
                x = x0 + sqrt(spacing^2 / (m(l)^2 + 1));
            else
                x = x0 - sqrt(spacing^2 / (m(l)^2 + 1));
            end
            y = m(l) * x - m(l) * x0 + y0;
            x0 = x;
            y0 = y;
            new_pts(k, :) = [x, y];
            k = k + 1;
        end
    end
    l = l + 1;
end

if(isempty(new_pts))
    pv = pv(1:end-1, :);
else
    pv = [pv_orig; new_pts];
end

end

