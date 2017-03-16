% Q3, HW 4
clear all
p_max = 1.0; u_max = 1.0;
p_left = p_max/2; p_right = 0.0;
dx = 4/400;
dt = 0.8 * dx / u_max;
XL = -2; X = 4; T = 28;

mesh = XL:dx:X;
time = 0:dt:T;

% number of time steps while red, equal to number of time steps while green
period = ceil(1/dt);
red1 = false; t_switch1 = 0;
red2 = true;

backsize = length(mesh(mesh < 0.0));
frontsize = length(mesh) - backsize;
back2 = length(mesh(mesh < 0.15));
front2 = length(mesh) - back2;

% bounding nodes for lights
light1 = [backsize, backsize + 1];
light2 = [back2, back2 + 1];

K = 0:9;
Loc = floor(back2 + 100);
averages = [];

for k = K
tau = floor(k * period / 10); % number of time steps to delay the second light
    
    pG = cell(length(time), 1);
    for i = 1:length(time)
        pG{i} = [0.0 .* ones(1, length(mesh))];
    end
    
    switch_flag = 0;
    for t = 1:length(time)        
        F_leftG = zeros(size(mesh)); F_rightG = zeros(size(mesh));
        
        switch_flag = switch_flag + 1;
        if mod(t, period) == 0
            switch_flag = 0;
            red1 = ~red1;
        end

        if switch_flag == tau
            red2 = ~red2;
        end

        for i = 1:(length(mesh) - 1)
            if i == 1
                left_pointG = p_left;
            else
                left_pointG = pG{t}(i - 1);
            end

            [F_leftG(i), F_rightG(i)] = Godunov2(left_pointG, pG{t}(i), pG{t}(i+1), p_max, u_max);

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
    
    pmatrix = [];
    for tt = 1:length(time)
        pmatrix(tt, :) = pG{tt};
    end
    
    % compute average flow
    averages = [averages, sum(flux(pmatrix(25*period:27*period, Loc), p_max, u_max))./ (2*period)];
end

plot(K, averages)
xlabel('k')
ylabel('Average Capacity')
ylim([0.08, 0.13])