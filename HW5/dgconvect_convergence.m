%function [errors, slopes] = dgconvect_convergence

T = 1;                      % end simulation time
dt = 2e-4;                  % time step size
P = [1, 2, 4, 8, 16];       % polynomial orders for convergence study
N = [16, 32, 48, 64, 80, 96, 112, 128, 256]; % numbers of nodes

errors = zeros(length(P), length(P));

for p = 1:length(P)
    en = 1;
    for n = N ./ P(p)
        [u, errors(p, en)] = dgconvect(n, P(p), T, dt);
        en = en + 1;
    end
end


slopes = [];
leg = cell(1, length(P)); leg_text = ''; % to hold the legend labels

for p = 1:length(P)
    % find the rates of convergence
    fit = polyfit(log(N ./ p), log(errors(p, :)), 1);
    slopes(p) = fit(1);
    loglog(N ./ p,errors(p, :), '*-')
    leg{p} = sprintf('p = %i, slope = %.2i', P(p), slopes(p));
    hold on
    
    if p == length(P)   
        leg_text = strcat(leg_text, sprintf(' leg{%i}', p));
    else
        leg_text = strcat(leg_text, sprintf(' leg{%i},', p));
    end
end

leg_text = strcat(strcat('legend(', leg_text), ')');

eval(leg_text)
xlabel('log(number of elements)')
ylabel('L2-norm')

%end