function [row] = row_number(p, p2)
% This function returns the row number of a point p in p2
tol = 1e-10;

for j = 1:length(p2(:,1))
    if (abs(p(1) - p2(j, 1)) < tol) && (abs(p(2) - p2(j, 2)) < tol)
        row = j;
    end
end
    
end

