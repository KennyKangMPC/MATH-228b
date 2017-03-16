% Q3, HW 4
clear all
p_max = 1.0; u_max = 1.0;
p_left = p_max/2; p_right = 0.0;
dx = 4/400;
dt = 0.8 * dx / u_max;
XL = -5; X = 5;
T = 30;

mesh = XL:dx:X;
time = 0:dt:T;

% number of time steps while red, equal to number of time steps while green
period = ceil(1/dt);
red1 = false;
red2 = true;

backsize = length(mesh(mesh < 0.0));
frontsize = length(mesh) - backsize;

back2 = length(mesh(mesh < 2.15));
front2 = length(mesh) - back2;

% bounding nodes for lights
light1 = [backsize, backsize + 1];
light2 = [back2, back2 + 1];

K = 9;
Loc = (backsize + 10):25:(backsize + 500);
averages = zeros(length(K), length(Loc));

r = 1;
for k = K
tau = floor(k * period / 10); % number of time steps to delay the second light

    pG = cell(length(time), 1);
    for i = 1:length(time)
        %pG{i} = [p_left .* ones(1, backsize), p_right .* ones(1, frontsize)];
        pG{i} = [p_left .* ones(1, back2), p_right .* ones(1, front2)];
    end

    for t = 1:length(time)
        F_leftG = zeros(size(mesh)); F_rightG = zeros(size(mesh));

        if mod(t, period) == 0
            red1 = ~red1;
        end

        if mod(t + tau, period) == 0
            red2 = ~red2;
        end

        for i = 1:(length(mesh) - 1)
            if i == 1
                left_pointG = p_left;
            else
                left_pointG = pG{t}(i - 1);
            end

            [F_leftG(i), F_rightG(i)] = Godunov(left_pointG, pG{t}(i), pG{t}(i+1), p_max, u_max);

            if (red1 && (i == light1(1)))
                F_rightG(i) = 0.0;
            end

            if (red1 && (i == light1(2)))
                F_leftG(i) = 0.0;
            end

            if (red2 && (i == light2(1)))
                F_rightG(i) = 0.0;
            end

            if (red2 && (i == light2(2)))
                F_leftG(i) = 0.0;
            end


        pG{t + 1}(i) = pG{t}(i) - (dt / dx) * (F_rightG(i) - F_leftG(i));
        end    
    end

    for t = 1:5:length(time)
        plot(mesh, pG{t})
        drawnow
    end
    
    % compute average flow
    l = 1;
    for location = Loc
        for t = (25*period):1:(27*period)
            averages(r, l) = averages(r, l) + flux(pG{t}(location), p_max, u_max)./period;
        end
        l = l + 1;
    end
    r = r + 1;
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
