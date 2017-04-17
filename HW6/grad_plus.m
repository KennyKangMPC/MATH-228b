function plus = grad_plus(phi, i, j, dx, dy, numx, numy)
    if i == 1
        a = max(phi(i, j)/dx, 0); % remove undefined variable
    else
        a = max((phi(i, j) - phi(i - 1, j))/dx, 0);
    end
    
    if i == numx
        b = min(-phi(i, j)/dx, 0); % remove undefined
    else
        b = min((phi(i + 1, j) - phi(i, j))/dx, 0); 
    end
    
    if j == 1
        c = max((phi(i, j)) / dy, 0); % remove undefined
    else
        c = max((phi(i, j) - phi(i, j - 1)) / dy, 0);
    end
    
    if j == numy
        d = min(-phi(i, j) / dy, 0);
    else
        d = min((phi(i, j + 1) - phi(i, j)) / dy, 0);
    end
        
    plus = sqrt(a.^2 + b.^2 + c.^2 + d.^2);
end