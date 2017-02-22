% solve for the shape functions
pts = [0,0; 0.5,0; 1,0; 0.5,0.5; 0,1; 0,0.5];

A = zeros(6, 6);
for i = 1:length(A(1,:))
    A(i, :) = [1, pts(i, 1), pts(i, 2), pts(i, 1).^2, pts(i, 2).^2, pts(i, 1) * pts(i, 2)];
end

for i = 1:length(A(1,:))
    b = zeros(6, 1);
    b(i) = 1;
    x = A\b;
    sprintf('Shape function %i: %f + %f*xi  +  %f*eta  +  %f*xi^2  +  %f*eta^2  +  %f*xi*eta', i, x(1), x(2), x(3), x(4), x(5), x(6))
end
