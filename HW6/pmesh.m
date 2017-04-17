function [p,t,e, data] = pmesh(pv, hmax, nref)
% PMESH  Delaunay refinement mesh generator.
% UC Berkeley Math 228B, Per-Olof Persson <persson@berkeley.edu>

data = struct('p', {}, 't', {}, 'e', {}, 'T', {}, 'R', {}, ...
              'A', {}, 'b', {}, 'u', {});

p = [];
for i = 1:size(pv,1)-1
  pp = pv(i:i+1,:);
  L = sqrt(sum(diff(pp,[],1).^2, 2));
  if L>hmax
    n = ceil(L / hmax);
    pp = interp1([0,1], pp, (0:n) / n);
  end
  p = [p; pp(1:end-1,:)];
end

while 1
  t = delaunayn(p);
  t = removeoutsidetris(p, t, pv);
  %tplot(p,t)

  area = triarea(p,t);
  [maxarea, ix] = max(area);
  if maxarea < hmax^2 / 2, break; end
  pc = circumcenter(p(t(ix,:),:));
  p(end+1,:) = pc;
end

e = boundary_nodes(t);

% save the first field
data(1).p = p;
data(1).t = t;
data(1).e = e;

for iref = 1:nref
    [pmid, ia] = edgemidpoints(p,t);
    p = [p; pmid];
    t = delaunayn(p);
    t = removeoutsidetris(p, t, pv);
    e = boundary_nodes(t);
    %tplot(p,t)

    % save all subsequent fields
    data(iref + 1).p = p;
    data(iref + 1).t = t;
    data(iref + 1).e = e;

    % compute the interpolation matrix to upgrade refinement iref
    old_pts = length(data(iref).p);
    new_pts = length(data(iref + 1).p);
    data(iref).T  = zeros(new_pts, old_pts);

    data(iref).T(1:old_pts, 1:old_pts) = eye(old_pts);
    num_tri = length(data(iref).t(:, 1)); % number of triangles from previous

    for i = 1:(new_pts - old_pts) % for each row in T that is _new_
        row = ia(i); % row from _entire_ pmid

        if row <= num_tri % on side 1-2
            subrow = row;
            data(iref).T(i + old_pts, data(iref).t(subrow, 1)) = 0.5;
            data(iref).T(i + old_pts, data(iref).t(subrow, 2)) = 0.5;
        elseif ((num_tri < row) && (row <= 2*num_tri)) % on side 2-3
            subrow = row - num_tri;
            data(iref).T(i + old_pts, data(iref).t(subrow, 2)) = 0.5;
            data(iref).T(i + old_pts, data(iref).t(subrow, 3)) = 0.5;
        else % on side 3-1
            subrow = row - 2*num_tri;
            data(iref).T(i + old_pts, data(iref).t(subrow, 3)) = 0.5;
            data(iref).T(i + old_pts, data(iref).t(subrow, 1)) = 0.5;
        end
    end
    
    % compute the reduction matrix as the transpose of T, with rows 
    % scaled to unity
    data(iref).R = transpose(data(iref).T);
    for i = 1:length(data(iref).R(:, 1)) % loop over the rows
        data(iref).R(i, :) = data(iref).R(i, :) ./ sum(data(iref).R(i, :));
    end   
end

function [pmid, ia] = edgemidpoints(p, t)
  
pmid = [(p(t(:,1),:) + p(t(:,2),:)) / 2;
        (p(t(:,2),:) + p(t(:,3),:)) / 2;
        (p(t(:,3),:) + p(t(:,1),:)) / 2];

[pmid, ia, ic] = unique(pmid, 'rows');


function a = triarea(p, t)

d12 = p(t(:,2),:) - p(t(:,1),:);
d13 = p(t(:,3),:) - p(t(:,1),:);
a = abs(d12(:,1).*d13(:,2) - d12(:,2).*d13(:,1)) / 2;

function t = removeoutsidetris(p, t, pv)

pmid = (p(t(:,1),:) + p(t(:,2),:) + p(t(:,3),:)) / 3;
isinside = inpolygon(pmid(:,1), pmid(:,2), pv(:,1), pv(:,2));
t = t(isinside,:);

function pc = circumcenter(p)

dp1 = p(2,:) - p(1,:);
dp2 = p(3,:) - p(1,:);

mid1 = (p(1,:) + p(2,:)) / 2;
mid2 = (p(1,:) + p(3,:)) / 2;
  
s = [-dp1(2),dp2(2); dp1(1),-dp2(1)] \ [-mid1 + mid2]';
pc = mid1 + s(1) * [-dp1(2),dp1(1)];