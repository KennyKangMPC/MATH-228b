% mesh domain
pv = [0,0; 5,0; 5,1; 3.1,1; 3.1,0.2; 2.9,0.2; 2.9,1; 2.1,1; 2.1,0.2; 1.9,0.2; 1.9,1.0; 0,1; 0,0];
[p, t, e] = pmesh(pv, 0.2, 2);
%tplot(p, t)

% find the mesh edges
[ein, eout, ewall] = waveguide_edges(p, t);

H = [];
Wave = 6:0.01:6.5;

% solve femhelmholtz
for wave = Wave

    [K, M, Bin, Bout, bin] = femhelmholtz(p, t, ein, eout, wave);
    Kk = K - (wave.^2) .* M + 1i .* wave .* (Bin + Bout);
    Ff = bin .* 2 .* 1i .* wave;
    a = Kk\Ff;

    % compute complex conjugate transpose
    a_h = transpose(conj(a));

    % compute H(u)
    H = [H, a_h * Bout * a];
end

semilogy(Wave, real(H), 'o-')
xlabel('k')
ylabel('log(H(u))')
grid on
%tplot(p, t, real(a))