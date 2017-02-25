% script for convergence study of femhelmholtz
clear all

wave = 6;
imag = 1i;
hmax = 0.3;
nrefmax = 4;
errors = zeros(1, nrefmax);
tol = 1e-10;

pv = [0,0; 5,0; 5,1; 0,1; 0,0];

for nref = 0:(nrefmax - 1)
    [p, t, e] = pmesh(pv, hmax, nref);
    [ein, eout, ewall] = waveguide_edges(p, t);
    [K, M, Bin, Bout, bin] = femhelmholtz(p, t, ein, eout);
    Kk = K - (wave.^2) .* M + imag .* wave .* (Bin + Bout);
    Ff = bin .* 2 .* imag .* wave;
    a = Kk\Ff;
    
    u_exact = cos(wave .* p(:,1));

    tplot(p, t, real(a) - u_exact)
    k = 1;
    error = [];

    % loop over all the points in the ref soln to find the maximum error
    for i = 1:length(p)
        error(k) = abs(real(a(i)) - u_exact(i));
        k = k + 1;
    end

    % determine the max-norm (maximum element in error())
    errors(nref + 1) = max(error);
end

sizes = 1./hmax .^ (0:(nrefmax - 1));
loglog(sizes, errors, 'o-')
xlabel('1/hmax', 'FontSize', 16)
ylabel('Max-Norm of the Error', 'FontSize', 16)
grid on

for i = 1:(length(errors) - 1)
    rate(i) = log2(errors(i)) - log2(errors(i+1));
end

saveas(gcf, 'error', 'png')
