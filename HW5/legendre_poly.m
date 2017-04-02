function [y, dy] = legendre_poly(x, p)
%LEGENDRE_POLY  Legendre polynomials and derivatives of degree p at nodes x

% UC Berkeley Math 228B, Per-Olof Persson <persson@berkeley.edu>
    
    z = 0*x;
    y = [z+1, x, z*ones(1,p-1)];
    dy = [z, z+1, z*ones(1,p-1)];
    for i = 1:p-1
        y(:,i+2) = ((2*i+1)*x.*y(:,i+1) - i*y(:,i)) / (i+1);
        dy(:,i+2) = ((2*i+1)*(x.*dy(:,i+1) + y(:,i+1)) - i*dy(:,i)) / (i+1);
    end
end
