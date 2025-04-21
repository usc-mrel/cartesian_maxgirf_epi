function b = siemens_B(r, theta, phi, R0, alpha, beta)

N = length(r); % number of voxels
nmax = size(alpha,1) - 1;
b = zeros(N, 1, 'double');
for n = 0:nmax
    f = R0 * (r / R0).^n;
    % Compute the associated Legendre functions of degree n and order m = 0, 1, ..., n
    p_unnorm = legendre(n, cos(theta)).'; % N x n+1
    for m = 0:n
        f2 = alpha(n + 1, m + 1) * cos(m * phi) + beta(n + 1, m + 1) * sin(m * phi);
        ptemp = p_unnorm(:,m + 1);
        normfact = 1;
        if m > 0
            normfact = (-1)^m * sqrt((2 * n + 1) * factorial(n - m) / (2 * factorial(n + m)));
        end
        p = normfact * ptemp;
        b = b + f .* p .* f2;
    end
end
end
