% Question 2, MATH-228b HW 1
clear all

xi = linspace(0, 1, 41);
eta = linspace(0, 1, 41);
T = 0.5;

[XI, ETA] = meshgrid(xi, eta);

a1 = @(xi) 2*xi^3 - 3*xi^2 + 1;
a2 = @(xi) 3*xi^2 -2*xi^3;
a3 = @(xi) xi^3 - 2*xi^2 + 1;
a4 = @(xi) xi^3 - xi^2;

b1 = @(eta) 2*eta^3 - 3*eta^2 + 1;
b2 = @(eta) 3*eta^2 -2*eta^3;
b3 = @(eta) eta^3 - 2*eta^2 + 1;
b4 = @(eta) eta^3 - eta^2;

leftx = @(xi,eta) 0;
lefty = @(xi,eta) eta;
rightx = @(xi,eta) 1;
righty = @(xi,eta) 1.4*eta;
bottomx = @(xi,eta) xi;
bottomy = @(xi,eta) (1-cos(2*pi*xi))/5;
topx = @(xi,eta) xi;
topy = @(xi,eta) 1 + (1-cos(pi*xi))/5;

k = 1;
for dxi = linspace(0, 1, 41)
    for deta = linspace(0, 1, 41)
        x(k) = a1(dxi)*leftx(dxi, deta) + a2(dxi)*rightx(dxi, deta) +...
            b1(deta)*bottomx(dxi, deta) + b2(deta)*topx(dxi, deta) - ...
            (a1(dxi)*b1(deta)*leftx(0, 0) + a1(dxi)*b2(deta)*leftx(0, 1) ...
            + a2(dxi)*b1(deta)*rightx(1,0) + a2(dxi)*b2(deta)*rightx(1,1));
        y(k) = a1(dxi)*lefty(dxi, deta) + a2(dxi)*righty(dxi, deta) +...
            b1(deta)*bottomy(dxi, deta) + b2(deta)*topy(dxi, deta) - ...
            (a1(dxi)*b1(deta)*lefty(0, 0) + a1(dxi)*b2(deta)*lefty(0, 1) ...
            + a2(dxi)*b1(deta)*righty(1,0) + a2(dxi)*b2(deta)*righty(1,1));
        k = k + 1;
    end
end

%x = unique(x);
%y = unique(y);
plot(x, y, 'o')




