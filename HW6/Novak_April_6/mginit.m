function data = mginit(pv, hmax, nref)

[p, t, e, data] = pmesh(pv, hmax, nref);

% determine the A and b matrices at the finest mesh by calling fempoi
[data(nref + 1).u, data(nref + 1).A, data(nref + 1).b] = ...
    fempoi(data(nref + 1).p, data(nref + 1).t, data(nref + 1).e);

% reduce the finest-mesh A and b 
for i = fliplr(1:nref)
    data(i).A = data(i).R * data(i + 1).A * data(i).T;
    data(i).b = data(i).R * data(i + 1).b;
end

end

