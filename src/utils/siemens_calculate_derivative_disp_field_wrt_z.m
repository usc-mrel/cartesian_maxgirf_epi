function dDeltak_dz = siemens_calculate_derivative_disp_field_wrt_z(r, theta, phi, R0, alpha, beta, Gref)
% Calculate the Cartesian derivative of a displacement field with respect to z

N = length(r); % number of voxels
dDeltak_dz = zeros(N, 1, 'single');

mpsi = atan2(beta, alpha);
nmax = size(alpha,1) - 1;
sin_theta = sin(theta);
cos_theta = cos(theta);

for n = 0:nmax
    f = (1 / R0)^n;
    r_n_minus_1 = r.^(n - 1);

    %----------------------------------------------------------------------
    % Compute the unnormalized associated Legendre functions of
    % order n and degree m = 0, 1, ..., n
    %----------------------------------------------------------------------
    p_unnorm = legendre(n, cos_theta).'; % N x n+1
    p_unnorm = bsxfun(@times, p_unnorm, (-1).^(0:n));

    for m = 0:n
        %------------------------------------------------------------------
        % Calculate the normalization factor
        %------------------------------------------------------------------
        N_nm = 1;
        if m > 0
            N_nm = (-1)^m * sqrt((2 * n + 1) * factorial(n - m) / (2 * factorial(n + m)));
        end

        %------------------------------------------------------------------
        % Calculate the unnormalized associated Legendre functions
        %------------------------------------------------------------------
        p_nm = p_unnorm(:,m+1);

        if m + 1 > n
            p_nm_plus_1 = zeros(N, 1, 'single');
        else
            p_nm_plus_1 = p_unnorm(:,m+1+1);
        end

        %------------------------------------------------------------------
        % Calculate d T_{n,m} / dz
        %------------------------------------------------------------------
        dTnm_dz = r_n_minus_1 .* cos(m * phi - mpsi(n+1,m+1)) .* ((n - m) * cos_theta .* p_nm + sin_theta .* p_nm_plus_1);

        %------------------------------------------------------------------
        % Calculate the summation
        %------------------------------------------------------------------
        dDeltak_dz = dDeltak_dz + (f * sqrt(alpha(n+1,m+1)^2 + beta(n+1,m+1)^2)) * N_nm * (-1)^m * dTnm_dz;
    end
end

dDeltak_dz = dDeltak_dz / Gref;

end
