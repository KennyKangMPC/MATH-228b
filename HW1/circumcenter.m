function [pt, max_area] = circumcenter(T, pv)
% finds the circumcenter of the largest triangle
% compute triangle areas using Heron's formula
for i = 1:length(T(:,1))
    A = [pv(T(i, 1), 1), pv(T(i, 1), 2)];
    B = [pv(T(i, 2), 1), pv(T(i, 2), 2)];
    C = [pv(T(i, 3), 1), pv(T(i, 3), 2)];
    area(i) = (A(1) * (B(2) - C(2)) + B(1) * (C(2) - A(2)) + C(1) * (A(2) - B(2))) / 2;
end

% find the largest triangle, and store the number in refine
max_area = area(1);
refine = 1;
for i = 1:length(area)
    if area(i) > max_area
        max_area = area(i);
        refine = i;
    end
end

sprintf('Refining triangle %i', refine)

A = [pv(T(refine, 1), 1), pv(T(refine, 1), 2)];
B = [pv(T(refine, 2), 1), pv(T(refine, 2), 2)];
C = [pv(T(refine, 3), 1), pv(T(refine, 3), 2)];

AB = [(A(1) + B(1)) / 2, (A(2) + B(2)) / 2];
m_AB = (B(2) - A(2)) / (B(1) - A(1));

AC = [(A(1) + C(1)) / 2, (A(2) + C(2)) / 2];
m_AC = (C(2) - A(2)) / (C(1) - A(1));

BC = [(B(1) + C(1)) / 2, (B(2) + C(2)) / 2];
m_BC = (C(2) - B(2)) / (C(1) - B(1));

if m_BC == 0
    m_AB = -1/m_AB;
    m_AC = -1/m_AC;
    mat = [1, -m_AC; 1, -m_AB];
    b = [-m_AC * AC(1) + AC(2); -m_AB * AB(1) + AB(2)];
elseif m_AC == 0
    m_AB = -1/m_AB;
    m_BC = -1/m_BC;
    mat = [1, -m_BC; 1, -m_AB];
    b = [-m_BC * BC(1) + BC(2); -m_AB * AB(1) + AB(2)];
else
    m_AC = -1/m_AC;
    m_BC = -1/m_BC;
    mat = [1, -m_BC; 1, -m_AC];
    b = [-m_BC * BC(1) + BC(2); -m_AC * AC(1) + AC(2)];
end

circumcenter = transpose(mat\b);
pt(1) = circumcenter(2);
pt(2) = circumcenter(1);

end

