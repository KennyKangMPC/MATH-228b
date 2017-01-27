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
for i = 1:1
    if pv(i + 1, 1) == pv(i, 1)
        % vertical line
    else
        % not a vertical line
        x0 = pv(i, 1);
        y0 = pv(i, 2);
        m = (pv(i + 1, 2) - y0) / (pv(i + 1, 1) - x0);
        % determine number of points to place to get close to hmax
        dist = sqrt((pv(i + 1, 2) - pv(i, 2))^2 + (pv(i + 1, 1) - pv(i, 1))^2);
        num_divs = ceil(dist/hmax);
        spacing = dist/num_divs;
        for j = 1:(num_divs - 1)
           x = sqrt(spacing^2 / (m^2 + 1)) + x0;
           y = m * x - m * x0 + y0;
           x0 = x;
           y0 = y;
           new_pts(j, :) = [x, y];
        end
    end
end
