function [T_outside_removed] = delete_outside_triangles(T, pv, pv_orig)
% delete triangles outside the domain

s = 1;
for i = 1:length(T(:,1))
    A = [pv(T(i, 1), 1), pv(T(i, 1), 2)];
    B = [pv(T(i, 2), 1), pv(T(i, 2), 2)];
    C = [pv(T(i, 3), 1), pv(T(i, 3), 2)];
    ptA = [(A(1) + B(1)) / 2, (A(2) + B(2)) / 2];
    ptB = [(A(1) + C(1)) / 2, (A(2) + C(2)) / 2];
    ptC = [(C(1) + B(1)) / 2, (C(2) + B(2)) / 2];
    A_in = inpolygon(ptA(1), ptA(2), pv_orig(:,1), pv_orig(:,2));
    B_in = inpolygon(ptB(1), ptB(2), pv_orig(:,1), pv_orig(:,2));
    C_in = inpolygon(ptC(1), ptC(2), pv_orig(:,1), pv_orig(:,2));

    if A_in && B_in && C_in
        T_outside_removed(s, :) = T(i, :);
        s = s + 1;
    end
end


end

