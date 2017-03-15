% Q1, HW 4
clear all
p_max = 1.0;
u_max = 1.0;
p_left = 0.8;
p_right = 0.0;
dx = 4/400;
dt = 0.8 * dx / u_max;
XL = -5;
X = 5;
T = 2;

mesh = XL:dx:X;
time = 0:dt:T;

backsize = length(mesh(mesh < 0.0));
frontsize = length(mesh) - backsize;

pG = cell(length(time), 1);
pR = cell(length(time), 1);
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
        
        flux_leftG = flux(left_pointG);
        flux_midG = flux(pG{t}(i));
        flux_rightG = flux(pG{t}(i + 1));
        
        flux_leftR = flux(left_pointR);
        flux_midR = flux(pR{t}(i));
        flux_rightR = flux(pR{t}(i + 1));
        
        [F_leftG(i), F_rightG(i)] = Godunov(flux_leftG, flux_midG, flux_rightG, left_pointG, pG{t}(i), pG{t}(i+1));
        [F_leftR(i), F_rightR(i)] = Roe(flux_leftR, flux_midR, flux_rightR, left_pointR, pR{t}(i), pR{t}(i+1));
        
    pG{t + 1}(i) = pG{t}(i) - (dt / dx) * (F_rightG(i) - F_leftG(i));
    pR{t + 1}(i) = pR{t}(i) - (dt / dx) * (F_rightR(i) - F_leftR(i));
    end    
end

for t = 1:length(time)
   plot(mesh, pG{t}, 'g', mesh, pR{t}, 'r');
   ylim([0, 0.8])
   drawnow
end
