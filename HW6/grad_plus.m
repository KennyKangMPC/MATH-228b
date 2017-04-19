function plus = grad_plus(phi, i, j, dx, dy, numx, numy)
    big = 1000;
    if i == 1
        a = max(phi(i, j)/dx - big, 0) ^ 2;
    else
        a = max((phi(i, j) - phi(i - 1, j))/dx, 0) ^ 2;
    end
    
    if i == numx
        b = min(big - phi(i, j)/dx, 0) ^ 2; % remove undefined
    else
        b = min((phi(i + 1, j) - phi(i, j))/dx, 0) ^ 2; 
    end
    
    if j == 1
        c = max((phi(i, j) - big) / dy, 0) ^ 2; % remove undefined
    else
        c = max((phi(i, j) - phi(i, j - 1)) / dy, 0) ^ 2;
    end
    
    if j == numy
        d = min(big -phi(i, j) / dy, 0) ^ 2;
    else
        d = min((phi(i, j + 1) - phi(i, j)) / dy, 0) ^ 2;
    end
        
    plus = sqrt(a + b + c + d);
end