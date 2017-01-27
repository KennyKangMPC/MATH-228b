% MATH 228b - HW1, question 4
clear all

hmax = 0.2;
% number of uniform refinements
nref = 0;
pv = [0,0; 1,0; .5,.5; 1,1; 0,1; 0,0];

% put nodes along the original edges of the polygon
%plot(pv(:,1), pv(:,2))

% find line slope
%for i = 1:(size(pv(:,1)) - 1)
k = 1;
for i = 1:5
    x0 = pv(i, 1);
    x_next = pv(i + 1, 1);
    y0 = pv(i, 2);
    y_next = pv(i + 1, 2);
    
    if x_next == x0
        % vertical line
        disp('vertical line')
        dist = abs(y_next - y0);
        num_divs = ceil(dist/hmax)
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
    else
        % not a vertical line
        m = (y_next - y0) / (x_next - x0);
        
        % determine number of points to place to get close to hmax
        dist = sqrt((y_next - y0)^2 + (x_next - x0)^2);
        num_divs = ceil(dist/hmax);
        spacing = dist/num_divs;
        
        % compute the new points
        for j = 1:(num_divs - 1)
            if x_next > x0
                x = x0 + sqrt(spacing^2 / (m^2 + 1));
            else
                x = x0 - sqrt(spacing^2 / (m^2 + 1));
            end
                y = m * x - m * x0 + y0;
                x0 = x;
                y0 = y;
                new_pts(k, :) = [x, y];
                k = k + 1;
        end
    end
end
plot(pv(:,1), pv(:,2), 'o')
hold on
plot(new_pts(:,1), new_pts(:,2), 'x')