function divF = spectral_divergence(Fx, Fy, h)

N = length(Fx(1, :));
L = N * h;
k = ones(N, 1) * [0:N/2-1 0 -N/2+1:-1];

% compute x derivative
dx = fft(Fx, [], 2);
w_hat_x = 1i * k .* dx;
wx = real(ifft(w_hat_x, [], 2));

% compute y derivative
dy = fft(Fy, [], 1);
w_hat_y = 1i * k' .* dy;
wy = real(ifft(w_hat_y, [], 1));

divF = (wx + wy) * L / (2*pi);
