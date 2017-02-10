function [pv] = delete_outside(pv, pv_orig)
% delete points outside the domain

j = 1;

for i = 1:size(pv(:,1))
    in = inpolygon(pv(i,1), pv(i,2), pv_orig(:,1), pv_orig(:,2));
    
    if in == 0        % generic point outside the domain
    else              % inside the domain
        pv_inside(j,:) = pv(i, :);
        j = j + 1;
    end
end

pv = pv_inside;

end

