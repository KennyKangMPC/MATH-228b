function [x, w] = gauss_quad(p)
%GAUSS_QUAD  Gaussian quadrature on [-1,1] for given degree of precision

% UC Berkeley Math 228B, Per-Olof Persson <persson@berkeley.edu>
    
    n = ceil((p+1)/2);
    b = 1:n-1;
    b = b ./ sqrt(4*b.^2-1);
    [Q,D] = eig(diag(b,1) + diag(b,-1));
    x = diag(D);
    w = 2*(Q(1,:).^2)';
end
