% mginit
clear all

% specify the mesh
pv = [0,0; 1,0; 1,1; 0,1; 0,0];
[p, t, e, data, pmid, ia] = pmesh(pv, 0.75, 1);

ref = 1;

old_pts = length(data(ref).p);
new_pts = length(data(ref + 1).p);
interp  = zeros(new_pts, old_pts);

interp(1:old_pts, 1:old_pts) = eye(old_pts);
num_tri = length(data(ref).t(:, 1)); % number of triangles from previous

for i = 1:(new_pts - old_pts) % for each row in interp that is _new_
    row = ia(i); % row from _entire_ pmid
    
    if row <= num_tri % on side 1-2
        subrow = row;
        interp(i + old_pts, data(ref).t(subrow, 1)) = 0.5;
        interp(i + old_pts, data(ref).t(subrow, 2)) = 0.5;
    elseif ((num_tri < row) && (row <= 2*num_tri)) % on side 2-3
        subrow = row - num_tri;
        interp(i + old_pts, data(ref).t(subrow, 2)) = 0.5;
        interp(i + old_pts, data(ref).t(subrow, 3)) = 0.5;
    else % on side 3-1
        subrow = row - 2*num_tri;
        interp(i + old_pts, data(ref).t(subrow, 3)) = 0.5;
        interp(i + old_pts, data(ref).t(subrow, 1)) = 0.5;
    end
end

% test out the interp matrix acting on a linear function z = x + y
x = data(ref).p(:, 1);
y = data(ref).p(:, 2);g
funct = x.*x; 

%tplot(data(ref).p, data(ref).t, funct)

funct_interp = interp * funct;
tplot(data(ref + 1).p, data(ref + 1).t, funct_interp)
% perform the solve to get the true solution u0
%[u_true, A, b] = fempoi(p, t, e);
%tplot(p, t, u_true)

%[u, res] = gauss_seidel(A, b, 0 .* u_true, 100);
%semilogy(0:niter, res)

