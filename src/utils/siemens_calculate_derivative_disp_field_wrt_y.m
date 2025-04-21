function dDeltak_dy = siemens_calculate_derivative_disp_field_wrt_y(r, theta, phi, R0, alpha, beta, Gref)
% Calculate the Cartesian derivative of a displacement field with respect to y

N = length(r); % number of voxels
dDeltak_dy = zeros(N, 1, 'single');

mpsi = atan2(beta, alpha);
nmax = size(alpha,1) - 1;
cos_theta = cos(theta);

for n = 0:nmax
    f = (1 / R0)^n;
    r_n_minus_1 = r.^(n - 1);

    %----------------------------------------------------------------------
    % Compute the unnormalized associated Legendre functions of
    % order n and degree m = 0, 1, ..., n
    %----------------------------------------------------------------------
    n_prime = n - 1;
    if n_prime < 0
        %fprintf('n = %2d: n_prime = n - 1 = %2d < 0: n_prime = -n_prime - 1 = %2d\n', n, n_prime, -n_prime - 1);
        n_prime = -n_prime - 1;
    else
        %fprintf('n = %2d: n_prime = n - 1 = %2d >= 0:\n', n, n_prime);
    end
    p_unnorm = legendre(n_prime, cos_theta).'; % N x n+1
    p_unnorm = bsxfun(@times, p_unnorm, (-1).^(0:n_prime));

    for m = 0:n
        %------------------------------------------------------------------
        % Calculate the normalization factor
        %------------------------------------------------------------------
        N_nm = 1;
        if m > 0
            N_nm = (-1)^m * sqrt((2 * n + 1) * factorial(n - m) / (2 * factorial(n + m)));
        end

        %------------------------------------------------------------------
        % Calculate the Kronecker delta function (delta_{m,0})
        %------------------------------------------------------------------
        if m == 0
            delta = 1;
        else
            delta = 0;
        end

        %------------------------------------------------------------------
        % Calculate P_{n_prime,m_prime1}(cos(theta)) = P_{n-1,m+1}(cos(theta))
        %------------------------------------------------------------------
        m_prime1 = m + 1;
        if abs(m_prime1) > n_prime
            %fprintf('\tm = %2d: m_prime1 = m + 1 = %2d, abs(m_prime1)(%2d) > n_prime(%2d), none\n', m, m_prime1, abs(m_prime1), n_prime);
            p1 = zeros(N, 1, 'single');
        elseif m_prime1 < 0
            %fprintf('\tm = %2d: m_prime1 = m + 1 = %2d, m_prime1 < 0\n', m, m_prime1);
            m_prime1 = -m_prime1;
            p1 = (-1)^m_prime1 * factorial(n_prime - m_prime1) / factorial(n_prime + m_prime1) * p_unnorm(:,m_prime1+1);
        else
            %fprintf('\tm = %2d: m_prime1 = m + 1 = %2d, m_prime1(%2d) <= n_prime(%2d) \n', m, m_prime1, m_prime1, n_prime);
            p1 = p_unnorm(:,m_prime1+1);
        end

        %------------------------------------------------------------------
        % Calculate P_{n_prime,m_prime2}(cos(theta)) = P_{n-1,m-1}(cos(theta))
        %------------------------------------------------------------------
        m_prime2 = m - 1;
        if abs(m_prime2) > n_prime
            %fprintf('\tm = %2d: m_prime2 = m - 1 = %2d, abs(m_prime2)(%2d) > n_prime(%2d), none\n', m, m_prime2, abs(m_prime2), n_prime);
            p2 = zeros(N, 1, 'single');
        elseif m_prime2 < 0
            %fprintf('\tm = %2d: m_prime2 = m - 1 = %2d, m_prime2 < 0\n', m, m_prime2);
            m_prime2 = -m_prime2;
            p2 = (-1)^m_prime2 * factorial(n_prime - m_prime2) / factorial(n_prime + m_prime2) * p_unnorm(:,m_prime2+1);
        else
            %fprintf('\tm = %2d: m_prime2 = m - 1 = %2d, m_prime2(%2d) <= n_prime(%2d) \n', m, m_prime2, m_prime2, n_prime);
            p2 = p_unnorm(:,m_prime2+1);
        end

        %------------------------------------------------------------------
        % Calculate d T_{n,m} / dy
        %------------------------------------------------------------------
        dTnm_dy = (r_n_minus_1 / 2) .* (-(1 + delta) * p1 .* sin((m + 1) * phi - mpsi(n+1,m+1)) + ...
                                        -(1 - delta) * (n + m) * (n + m - 1) * p2 .* sin((m - 1) * phi - mpsi(n+1,m+1)));

        %------------------------------------------------------------------
        % Calculate the summation
        %------------------------------------------------------------------
        dDeltak_dy = dDeltak_dy + (f * sqrt(alpha(n+1,m+1)^2 + beta(n+1,m+1)^2)) * N_nm * (-1)^m * dTnm_dy;
    end
end

dDeltak_dy = dDeltak_dy / Gref;

end
