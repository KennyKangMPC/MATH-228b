% Q1, HW 4
clear all
p_max = 1.0; u_max = 1.0;
p_left = 0.8; p_right = 0.0;
dx = 4/400;
dt = 0.8 * dx / u_max;
XL = -5; X = 5;
T = 2;

mesh = XL:dx:X;
time = 0:dt:T;

backsize = length(mesh(mesh < 0.0));
frontsize = length(mesh) - backsize;

pG = cell(length(time), 1); pR = cell(length(time), 1);
for i = 1:length(time)
    pG{i} = [p_left .* ones(1, backsize), p_right .* ones(1, frontsize)];
    pR{i} = [p_left .* ones(1, backsize), p_right .* ones(1, frontsize)];
end

for t = 1:length(time)
    F_leftG = zeros(size(mesh)); F_rightG = zeros(size(mesh));
    F_leftR = zeros(size(mesh)); F_rightR = zeros(size(mesh));

    for i = 1:(length(mesh) - 1)
        if i == 1
            left_pointG = p_left;
            left_pointR = p_left;
        else
            left_pointG = pG{t}(i - 1);
            left_pointR = pR{t}(i - 1);
        end
        
        [F_leftG(i), F_rightG(i)] = Godunov(left_pointG, pG{t}(i), pG{t}(i+1), p_max, u_max);
        [F_leftR(i), F_rightR(i)] = Roe(left_pointR, pR{t}(i), pR{t}(i+1), p_max, u_max);
        
    pG{t + 1}(i) = pG{t}(i) - (dt / dx) * (F_rightG(i) - F_leftG(i));
    pR{t + 1}(i) = pR{t}(i) - (dt / dx) * (F_rightR(i) - F_leftR(i));
    end    
end

% plot the solutions
dc = 0.0;
Godunov = false; Roe = false; difference = false;
plottime = 1:50:length(time);
for t = plottime
    Gcolor = [1.0 - dc, 0.0, dc];
    if (Godunov)
        plot(mesh, pG{t}, 'Color', Gcolor)
        ylabel('Godunov Solution')
    elseif (Roe)
        plot(mesh, pR{t}, 'Color', Gcolor)
        ylabel('Roe Solution')
    elseif (difference)
        plot(mesh, pG{t} - pR{t}, 'Color', Gcolor)
        ylabel('Godunov solution - Roe solution')
    else
        plot(mesh, pG{t}, 'g', mesh, pR{t}, 'r')
        drawnow
    end
        
    hold on
    ylim([0, 0.8])
    xlim([-2, 2])
    xlabel('Spatial Domain')
    dc = dc + 0.1;
end
legend('t=0.0', 't=0.4', 't=0.8', 't=1.2', 't=1.6', 't=2.0')
