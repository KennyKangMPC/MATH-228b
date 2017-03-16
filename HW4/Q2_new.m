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
        
        if i == backsize && red
            F_rightG(i) = 0.0;
        end
        
        if i == backsize + 1 && red % just to the right of the light
            F_leftG(i) = 0.0;
        end
        
    pG{t + 1}(i) = pG{t}(i) - (dt / dx) * (F_rightG(i) - F_leftG(i));
    end    
end

for t = 1:5:(length(time) - 1/dt)
    plot(mesh, pG{t}, 'g')
    ylim([0, 1])
    drawnow
end
%legend('t=0.0', 't=0.4', 't=0.8', 't=1.2', 't=1.6', 't=2.0')
