% Question 2, MATH-228b HW 1
clear all

xi = linspace(0, 1, 41);
eta = linspace(0, 1, 41);
T = 0.5;
Topp = 0.;

[XI, ETA] = meshgrid(xi, eta);

a1 = 2*XI^3 - 3*XI^2 + 1;
a2 = 3*XI^2 -2*XI^3;
a3 = XI^3 - 2*XI^2 + 1;
a4 = XI^3 - XI^2;

b1 = 2*ETA^3 - 3*ETA^2 + 1;
b2 = 3*ETA^2 -2*ETA^3;
b3 = ETA^3 - 2*ETA^2 + 1;
b4 = ETA^3 - ETA^2;

x00 = [0, 0];
x01 = [0, 1];
x10 = [1, 0];
x11 = [1, 1.4]; % y-value defined separately

x_xi_0 = [XI, (1-cos(2*pi*XI))/5];
x_xi_1 = [XI, 1 + (1-cos(pi*XI))/5];
x_0_eta = [0*ETA, ETA];
x_1_eta = [0*ETA + 1, 1.4*ETA];
n_0_eta = [1, 0];
n_1_eta = [1, 0];
n_xi_0 = [0, 1];
n_xi_1 = [0, 1];

x = (a1*x_0_eta(1)) + (a2*x_1_eta(1)) + (a3*T*n_0_eta(1)) + (a4*T*n_1_eta(1));
x = x + (b1*x_xi_0(1)) + (b2*x_xi_1(1)) + (b3*T*n_xi_0(1)) + (b4*T*n_xi_1(1));
x = x - (a1*b1*x00(1)) - (a1*b2*x01(1)) - (a2*b1*x10(1)) - (a2*b2*x11(1));

x = x - (a1*b3*x00(1)*T*n_xi_0(1)) - (a1*b4*x01(1)*T*n_xi_1(1)) - (a2*b3*x10(1)*T*n_xi_0(1)) - (a2*b4*x11(1)*n_xi_1(1));
x = x - (b1*a3*x00(1)*T*n_0_eta(1)) - (b1*a4*x10(1)*T*n_1_eta(1)) - (b2*a3*x01(1)*T*n_0_eta(1)) - (b2*a4*x11(1)*T*n_1_eta(1));



y = (a1*x_0_eta(2)) + (a2*x_1_eta(2)) + (a3*T*n_0_eta(2)) + (a4*T*n_1_eta(2));                           % fixed eta
y = y + (b1*x_xi_0(2)) + (b2*x_xi_1(2)) + (b3*T*n_xi_0(2)) + (b4*T*n_xi_1(2));
y = y - (a1*b1*x00(2)) - (a1*b2*x01(2)) - (a2*b1*x10(2)) - (a2*b2*x11(2));

y = y - (a1*b3*x00(2)*T*n_xi_0(2)) - (a1*b4*x01(2)*T*n_xi_1(2)) - (a2*b3*x10(2)*T*n_xi_0(2)) - (a2*b4*x11(2)*n_xi_1(2));
y = y - (b1*a3*x00(2)*T*n_0_eta(2)) - (b1*a4*x10(2)*T*n_1_eta(2)) - (b2*a3*x01(2)*T*n_0_eta(2)) - (b2*a4*x11(2)*T*n_1_eta(2));




surf(x, y, zeros.*x)
