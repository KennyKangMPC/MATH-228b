clear all

%pv = [0,0; 2,0; 1.5,1; .5,1; 0,0];
pv = [0,0; 1,0; 1,1; 0,1; 0,0];
for iref = 1:3
    sprintf('Beginning refinement level %i', iref)
    data = mginit(pv, 0.5, iref);
    [u, res] = mgsolve(data, 2, 2, 1e-10, iref);
    semilogy(res, '*-'), hold on
end
legend('1', '2', '3')
hold off
