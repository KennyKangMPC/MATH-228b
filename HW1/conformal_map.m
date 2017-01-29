% Question 3, MATH 228b HW 1
x = linspace(0, 1, 21);
y = linspace(0, 2*pi, 81);

[X, Y] = meshgrid(x, y);

u = (6.*exp(2.*X)-13.*exp(X).*cos(Y)+6)./(9.*exp(2.*X)-12.*exp(X).*cos(Y)+4);
v = 5.*exp(X).*sin(Y)./(9.*exp(2.*X)-12.*exp(X).*cos(Y)+4);

surf(u, v, zeros.*u)
