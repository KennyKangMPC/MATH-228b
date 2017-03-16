% Q1, HW 4
clear all
p_max = 1.0; u_max = 1.0;
p_left = p_max/2; p_right = 0.0;
dx = 4/400;
dt = 0.8 * dx / u_max;
XL = -5; X = 5;
T = 50;

mesh = XL:dx:X;
time = 0:dt:T;

% number of time steps while red, equal to number of time steps while green
period = ceil(1/dt);
red = false;

backsize = length(mesh(mesh < 0.0));
frontsize = length(mesh) - backsize;

% bounding nodes for light 1
light1 = [backsize, backsize + 1];

pG = cell(length(time), 1);
for i = 1:length(time)
    pG{i} = [p_left .* ones(1, backsize), p_right .* ones(1, frontsize)];
end

for t = 1:length(time)
    F_leftG = zeros(size(mesh)); F_rightG = zeros(size(mesh));

    if mod(t, period) == 0
        red = ~red;
    end
    
    for i = 1:(length(mesh) - 1)
        if i == 1
            left_pointG = p_left;
        else
            left_pointG = pG{t}(i - 1);
        end
        
        [F_leftG(i), F_rightG(i)] = Godunov(left_pointG, pG{t}(i), pG{t}(i+1), p_max, u_max);
        
        if red
            if i == light1[i]
                F_rightG(i) = 0.0;
            elseif i == light1[2]
                F_leftG(i) = 0.0;
            else
            end
        end

    pG{t + 1}(i) = pG{t}(i) - (dt / dx) * (F_rightG(i) - F_leftG(i));
    end    
end

% plot the solutions
% plottime = floor([15.1, 15.5, 16, 16.5, 17] ./ dt);
% for t = plottime
%     plot(mesh, pG{t})
%     ylabel('Godunov Solution')
%     hold on
%     ylim([0, 1.0])
%     xlim([-2, 2])
%     xlabel('Spatial Domain')
% end
% legend('t=15.1', 't=15.5', 't=16.0', 't=16.5', 't=17.0')

% compute average flow
for location = (backsize + 100):10:(backsize + 500)
    sum = 0;
    for t = (15*period):1:(17*period)
        sum = sum + flux(pG{t}(location), p_max, u_max)./period;
    end
    sprintf('The estimated flow at location %d is %d', location*dx, sum)
end