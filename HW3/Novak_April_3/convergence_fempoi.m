% script for convergence study of fempoi2
clear all
pv = [0,0; 1,0; 1,1; 0,1; 0,0];
hmax = 0.3;
nrefmax = 4;

errors = zeros(1, nrefmax);

% create initial (very fine) mesh 
[p_ref, t, e] = pmesh(pv, hmax, nrefmax);
[p2, t2, e2] = p2mesh(p_ref, t);

% compute the reference solution
[a_ref] = fempoi2(p2, t2, e2);
tplot(p2, t2, a_ref)

for nref = 0:(nrefmax - 1)
    [p, t, e] = pmesh(pv, hmax, nref);
    [p, t, e] = p2mesh(p, t);
    [a] = fempoi2(p, t, e);
    
    if nref == 0
        P = length(p);
    end

    k = 1;
    error = max(abs(a_ref(1:P) - a(1:P)));

    % determine the max-norm (maximum element in error())
    errors(nref + 1) = error;
end

loglog(1./hmax .^ (0:(nrefmax - 1)), errors, 'o-')
xlabel('log(1/hmax)')
ylabel('Max-Norm of the Error')
grid on

for i = 1:(length(errors) - 1)
    rate(i) = log2(errors(i)) - log2(errors(i+1));
end