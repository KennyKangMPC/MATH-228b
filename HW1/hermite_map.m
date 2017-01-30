% Question 2, MATH-228b HW 1
clear all

xi = linspace(0, 1, 41);
eta = linspace(0, 1, 41);
T = 0.5;

[XI, ETA] = meshgrid(xi, eta);

a1 = 2*XI^3 - 3*XI^2+1;
a2 = 3*XI^2 -2*XI^3;
a3 = XI^3 - 2*XI^2 + 1;
a4 = XI^3 - XI^2;

b1 = 2*ETA^3 - 3*ETA^2+1;
b2 = 3*ETA^2 -2*ETA^3;
b3 = ETA^3 - 2*ETA^2 + 1;
b4 = ETA^3 - ETA^2;

%x_0_eta = [0; ETA];
%x_1_eta = [1; 1.4*ETA];
%dx_0_eta_dxi = 1;
%dx_1_eta_dxi = 1;

%x_xi_0 = [XI; (1 - cos(2*pi*XI))/5];
%x_xi_1 = [XI; 1 + (1 - cos(pi*XI))/5];
%dx_xi_0_deta = 1;
%dx_xi_1_deta = 1;

%x_0_0 = [0; 0];
%x_0_1 = [0; 1];
%x_1_0 = [0; 1];
%x_1_1 = [1; 1.4];

% conformal map without the boolean sum terms
x = (a1 * 0) + (a2 * 1) + (a3 * 1) + (a4 * 1);
x = x + (b1 * XI) + (b2 * XI) + (b3 * 1) + (b4 * 1);
%x = x - a2 * b2 * 1.0;

y = (a1 * ETA) + (a2 * 1.4 * ETA) + (a3 * 1) + (a4 * 1);
y = y + (b1 * (1 - cos(2*pi*XI))/5) + b2 * (1 + (1 - cos(pi*XI))/5) + (b3 * 1) + (b4 * 1);
%y = y - a1 * b2 * 1 - a2 * b1 * 1 - a2 * b2 * 1.4;



surf(x, y, zeros.*x)
