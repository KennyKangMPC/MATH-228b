%function [errors, slopes] = dgconvect_convergence
clear all
T = 1;                      % end simulation time
dt = 2e-4;                  % time step size
%P = [1, 2, 4, 8, 16];       % polynomial orders for convergence study
P = [1];
N = 16:16:256;              % numbers of nodes
k = 1e-3;                   % thermal conductivity

errors = zeros(length(P), length(P));

for p = 1:length(P)
    en = 1;
    for n = N ./ P(p)
        [u, errors(p, en)] = dgconvdiff(n, p, T, dt, k);
        en = en + 1;
    end
end


slopes = [];
leg = cell(1, length(P)); leg_text = ''; % to hold the legend labels

for p = 1:length(P)
    p_error = errors(p, :);
    % to remove points affected by rounding:
    %p_error = p_error(p_error > 1e-8);
    %N = N(1:length(p_error));
    % find the rates of convergence
    fit = polyfit(log(N), log(p_error), 1);
    slopes(p) = fit(1);
    loglog(N, p_error, '*-')
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
xlabel('log(number of nodes)')
ylabel('log(L2-norm)')

%end