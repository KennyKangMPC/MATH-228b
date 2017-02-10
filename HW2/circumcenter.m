function [pt] = circumcenter(A, B, C, pv)

AB = [(A(1) + B(1)) / 2, (A(2) + B(2)) / 2];
m_AB = (B(2) - A(2)) / (B(1) - A(1));

AC = [(A(1) + C(1)) / 2, (A(2) + C(2)) / 2];
m_AC = (C(2) - A(2)) / (C(1) - A(1));

BC = [(B(1) + C(1)) / 2, (B(2) + C(2)) / 2];
m_BC = (C(2) - B(2)) / (C(1) - B(1));

if abs(m_BC - 0) < 1e-3
    m_AB = -1/m_AB;
    m_AC = -1/m_AC;
    mat = [1, -m_AC; 1, -m_AB];
    b = [-m_AC * AC(1) + AC(2); -m_AB * AB(1) + AB(2)];
elseif abs(m_AC - 0) < 1e-3
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

