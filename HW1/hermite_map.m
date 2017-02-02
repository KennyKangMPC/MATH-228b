% Question 2, MATH-228b HW 1
clear all

T = 0.5;

a1 = @(xi) 2*xi^3 - 3*xi^2 + 1;
a2 = @(xi) 3*xi^2 -2*xi^3;
a3 = @(xi) xi^3 - 2*xi^2 + xi;
a4 = @(xi) xi^3 - xi^2;

b1 = @(eta) 2*eta^3 - 3*eta^2 + 1;
b2 = @(eta) 3*eta^2 -2*eta^3;
b3 = @(eta) eta^3 - 2*eta^2 + eta;
b4 = @(eta) eta^3 - eta^2;

leftx =     @(xi,eta) 0;
lefty =     @(xi,eta) eta;
rightx =    @(xi,eta) 1;
righty =    @(xi,eta) 1.4*eta;
bottomx =   @(xi,eta) xi;
bottomy =   @(xi,eta) (1-cos(2*pi*xi))/5;
topx =      @(xi,eta) xi;
topy =      @(xi,eta) 1 + (1-cos(pi*xi))/5;

leftnx =    @(xi,eta) 1;
leftny =    @(xi,eta) 0;
rightnx =   @(xi,eta) 1;
rightny =   @(xi,eta) 0;

bottomnx = @(xi,eta) (-2*pi*sin(2*pi*xi)/5)/sqrt((-2*pi*sin(2*pi*xi)/5)^2+1);
bottomny = @(xi,eta) (1)/sqrt((-2*pi*sin(2*pi*xi)/5)^2+1);
topnx = @(xi,eta) (-pi*sin(pi*xi)/5)/sqrt((-pi*sin(pi*xi)/5)^2+1);
topny = @(xi,eta) (1)/sqrt((-pi*sin(pi*xi)/5)^2+1);

k = 1;
for dxi = linspace(0, 1, 41)
    for deta = linspace(0, 1, 41)
        x(k) = a1(dxi)*leftx(dxi, deta) + a2(dxi)*rightx(dxi, deta) + b1(deta)*bottomx(dxi, deta) + b2(deta)*topx(dxi, deta) ... % original four terms
             + a3(dxi)*T*leftnx(dxi,deta) + a4(dxi)*T*rightnx(dxi,deta) + b3(deta)*T*bottomnx(dxi,deta) + b4(deta)*T*topnx(dxi,deta) ... % four derivative terms
              - (a1(dxi)*b1(deta)*leftx(0, 0) + a1(dxi)*b2(deta)*leftx(0, 1) + a2(dxi)*b1(deta)*rightx(1,0) + a2(dxi)*b2(deta)*rightx(1,1)) ... % subtract out corners
              - (a1(dxi)*b3(deta)*T*bottomnx(0,0)   + a1(dxi)*b4(deta)*T*topnx(0,1))  ...
              - (a2(dxi)*b3(deta)*T*bottomnx(1,0)   + a2(dxi)*b4(deta)*T*topnx(1,1))  ...
              - (b1(deta)*a3(dxi)*T*leftnx(0,0)     + b1(deta)*a4(dxi)*T*rightnx(1,0))  ...
              - (b2(deta)*a3(dxi)*T*leftnx(0,1)     + b2(deta)*a4(dxi)*T*rightnx(1,1));
        
        y(k) = a1(dxi)*lefty(dxi, deta) + a2(dxi)*righty(dxi, deta) + b1(deta)*bottomy(dxi, deta) + b2(deta)*topy(dxi, deta) ... % original four terms
             + a3(dxi)*T*leftny(dxi,deta) + a4(dxi)*T*rightny(dxi,deta) + b3(deta)*T*bottomny(dxi,deta) + b4(deta)*T*topny(dxi,deta) ... % four derivative terms
            - (a1(dxi)*b1(deta)*lefty(0, 0) + a1(dxi)*b2(deta)*lefty(0, 1) + a2(dxi)*b1(deta)*righty(1,0) + a2(dxi)*b2(deta)*righty(1,1)) ... % subtract out corners
            - (a1(dxi)*b3(deta)*T*bottomny(0,0)    + a1(dxi)*b4(deta)*T*topny(0,1)) ...
            - (a2(dxi)*b3(deta)*T*bottomny(1,0)    + a2(dxi)*b4(deta)*T*topny(1,1)) ...
            - (b1(deta)*a3(dxi)*T*leftny(0,0)      + b1(deta)*a4(dxi)*T*rightny(1,0)) ...
            - (b2(deta)*a3(dxi)*T*leftny(0,1)      + b2(deta)*a4(dxi)*T*rightny(1,1));
        k = k + 1;
    end
end

%x = unique(x);
%y = unique(y);
scatter(x, y, 'o')
saveas(gcf, 'hermite_plot', 'png')




