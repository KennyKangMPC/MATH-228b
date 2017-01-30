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

% conformal map without the boolean sum terms
x = (a1*0) + (a2*1) + (a3*T*1) + (a4*T*1);              % fixed eta
x = x + (b1*XI) + (b2*XI) + (b3*0) + (b4*0);            % fixed xi
x = x - (a1*b1*0) - (a1*b2*0) - (a2*b1*1) - (a2*b2*1);  % subtract corners

x = x - (a1*b3*0*0) - (a1*b4*0*0) - (a2*b3*1*0) - (a2*b4*1*0); % subtract mixed
x = x - (b1*a3*0*T) - (b1*a4*1*T) - (b2*a3*0*T) - (b2*a4*1*T); % subtract mixed

y = (a1*ETA) + (a2*1.4*ETA) + (a3*0) + (a4*0);                            % fixed eta
y = y + (b1*(1 - cos(2*pi*XI))/5) + (b2*(1 + (1 - cos(pi*XI))/5)) + (b3*T) + (b4*T); % fixed xi
y = y - (a1*b1*0) - (a1*b2*1) - (a2*b1*0) - (a2*b2*1.4);

y = y - (a1*b3*0*T) - (a1*b4*1*T) - (a2*b3*0*T) - (a2*b4*1.4*ETA*T);
y = y - (b1*a3*0*0) - (b1*a4*0*0) - (b2*a3*1*0) - (b2*a4*1.4*ETA*0);



surf(x, y, zeros.*x)
