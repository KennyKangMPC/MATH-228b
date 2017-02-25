function [LM] = changeLM(num_elem, LM, In, Out)
% change the LM for weird elements from 1 - 2 - 3 (2 - 3 is on boundary)
% to 2 - 3 - 1 (edge 2-3 is now the xe edge)

for elem = 1:num_elem
    % find the edges of the current element
    edges = [LM(elem,[1,2]); LM(elem,[2,3]); LM(elem,[3,1])];
    edges = sort(edges, 2);
    xe_edge = [LM(elem, 1), LM(elem, 2)];
    eta_edge = [LM(elem, 3), LM(elem, 1)];
    xe_edge = sort(xe_edge, 2);
    eta_edge = sort(eta_edge, 2);
    
    for edge = 1:length(edges(:,1))
        % is it on the In boundary?
        for in = 1:length(In(:, 1))
            if edges(edge, 1) == In(in, 1) && edges(edge, 2) == In(in, 2)                
                if edges(edge, 1) == xe_edge(1) && edges(edge, 2) == xe_edge(2)
                elseif edges(edge, 1) == eta_edge(1) && edges(edge, 2) == eta_edge(2)
                else
                    %sprintf('Changing LM for element %i', elem)
                    LM(elem, :) = [LM(elem, 2), LM(elem, 3), LM(elem, 1)];
                end
            end
        end
        
        % is it on the Out boundary?
        for in = 1:length(Out(:, 1))
            if edges(edge, 1) == Out(in, 1) && edges(edge, 2) == Out(in, 2)
                if edges(edge, 1) == xe_edge(1) && edges(edge, 2) == xe_edge(2)
                elseif edges(edge, 1) == eta_edge(1) && edges(edge, 2) == eta_edge(2)
                else
                    %sprintf('Changing LM for element %i', elem)
                    LM(elem, :) = [LM(elem, 2), LM(elem, 3), LM(elem, 1)];
                end
            end
        end
    end
end

end

