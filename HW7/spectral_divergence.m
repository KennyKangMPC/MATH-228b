function divF = spectral_divergence(Fx, Fy, h)
N = length(Fx);     % number of nodes
L = N * h;          % domain length

% compute x derivative
dx = fft(Fx);
w_hat_x = 1i * [0:N/2-1 0 -N/2+1:-1] .* dx;
wx = real(ifft(w_hat_x));

% compute y derivative
dy = fft(Fy);
w_hat_y = 1i * [0:N/2-1 0 -N/2+1:-1] .* dy;
wy = real(ifft(w_hat_y));

divF = wx + wy;