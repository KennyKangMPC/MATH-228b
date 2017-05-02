function w = spectral_filter(u)

N = length(u(1, :));
k0 = 3 * N / 8;
k = ones(N, 1) * [0:N/2-1 0 -N/2+1:-1];

% filter along the x direction
T = 1 ./ (1 + (k ./ k0) .^ 16);

u_hat = fft(u, [], 2);
w_hat = T .* u_hat;
w = real(ifft(w_hat, [], 2));

% then filter along the y direction using the x-filtered solution
T = 1 ./ (1 + (k' ./ k0) .^ 16);

u_hat = fft(w, [], 1);
w_hat = T .* u_hat;
w = real(ifft(w_hat, [], 1));