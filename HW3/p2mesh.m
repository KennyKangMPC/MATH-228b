% first, call pmesh to generate the starting mesh (works for linear
% elements)
pv = [0,0; 1,0; 1,1; 0,1; 0,0];
hmax = 0.15;
nref = 0.0;
[p,t,e] = pmesh(pv, hmax, nref);

% plot the starting mesh
tplot(p, t)
hold on

% add all of the midpoints
for i = 1:length(t(:,1))
    A = [p(t(i, 1), 1), p(t(i, 1), 2)];
    B = [p(t(i, 2), 1), p(t(i, 2), 2)];
    C = [p(t(i, 3), 1), p(t(i, 3), 2)];

    % add the three midpoints
    p = [p; (A(1) + B(1))/2, (A(2) + B(2))/2];
    p = [p; (A(1) + C(1))/2, (A(2) + C(2))/2];
    p = [p; (C(1) + B(1))/2, (C(2) + B(2))/2];
end

% p2 contains both the original mesh points and all of the midpoints
% (none duplicated)
p2 = unique(p, 'rows');
scatter(p2(:,1), p2(:,2), 'o')