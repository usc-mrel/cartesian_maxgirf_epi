function outp = cartesian_maxgirf_2d_normal(inp, sens, mask, p1, p2, iflag, eps, support_mask, U, V, L, initial_phase, lambda)

%% Start a stopwatch timer
start_time = tic;

%% Declare a persistent variable
persistent cg_iter;
if isempty(cg_iter)
    cg_iter = 0;
else
    cg_iter = cg_iter + 1;
end

%% Get imaging parameters
[Nkx,Nky,Nkz,Nc] = size(sens);

%--------------------------------------------------------------------------
% Calculate the number of k-space samples
%--------------------------------------------------------------------------
Nk = size(U,1);

%--------------------------------------------------------------------------
% Calculate the total number of voxels
%--------------------------------------------------------------------------
N = Nkx * Nky * Nkz;

%--------------------------------------------------------------------------
% Calculate the number of voxels within a support mask
%--------------------------------------------------------------------------
N_support = length(find(support_mask));

%% Calculate scale factors
type1_scale_factor = 1 / sqrt(Nkx * Nky * Nkz);
type2_scale_factor = 1 / sqrt(Nkx * Nky * Nkz);

%% Calculate P * m
Pm = bsxfun(@times, initial_phase, reshape(inp, [Nkx Nky Nkz])); % Nkx x Nky x Nkz

%% Calculate S * P * m
SPm = bsxfun(@times, sens, Pm); % Nkx x Nky x Nkz x Nc

%% Calculate d = sum_{ell=1}^{L} (I_{Nc} kron U_{ell}) * (I_{Nc} kron R_{Omega} * F) * (I_{Nc} kron V_{ell}^H) * S * P * m
d = complex(zeros(Nk, Nc, 'single'));

for ell = 1:L
    %----------------------------------------------------------------------
    % Calculate (I_{Nc} kron V_{ell}^H) * S * m
    %----------------------------------------------------------------------
    VhSPm = bsxfun(@times, conj(V(:,:,:,ell)), SPm); % Nkx x Nky x Nkz x Nc
    VhSPm = reshape(VhSPm, [N Nc]); % N x Nc

    %----------------------------------------------------------------------
    % Perform type-1 NUFFT (nonuniform, image space => uniform, k-space)
    % Siemens: k-space <=> image space
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s:(PCG=%d)(ell=%2d/%2d) Performing type-1 NUFFT... ', datetime, cg_iter, ell, L);
    cj = reshape(VhSPm((support_mask > 0),:), [N_support Nc]);
    FVhSPm = reshape(type1_scale_factor * finufft2d1(p1, p2, cj, -iflag, eps, Nkx, Nky), [N Nc]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Calculate (I_{Nc} kron R_{Omega} * F) * (I_{Nc} kron V_{ell}^H) * S * m
    %----------------------------------------------------------------------
    RFVhSPm = FVhSPm((mask > 0),:); % Nk x Nc

    %----------------------------------------------------------------------
    % Calculate (I_{Nc} kron U_{ell}) * (I_{Nc} kron R_{Omega} * F) * (I_{Nc} kron V_{ell}^H) * S * m
    %----------------------------------------------------------------------
    d = d + bsxfun(@times, U(:,ell), RFVhSPm); % Nk x Nc
end

%% Calculate sum_{ell=1}^{L} (I_{Nc} kron V_{ell}) * (I_{Nc} kron F^H * R_{Omega}^H) * (I_{Nc} kron U_{ell}^H) * d
VFhRhUhd = complex(zeros(Nkx, Nky, Nkz, Nc, 'single'));

for ell = 1:L
    %----------------------------------------------------------------------
    % Calculate (I_{Nc} kron U_{ell}^H) * d
    %----------------------------------------------------------------------
    Uhd = bsxfun(@times, conj(U(:,ell)), d); % Nk x Nc

    %----------------------------------------------------------------------
    % Calculate (I_{Nc} kron F^H * R_{Omega}^H) * (I_{Nc} kron U_{ell}^H) * d
    %----------------------------------------------------------------------
    RhUhd = complex(zeros(N, Nc, 'single'));
    RhUhd((mask > 0),:) = Uhd; % Nk x Nc => N x Nc
    RhUhd = reshape(RhUhd, [Nkx Nky Nc]);

    %----------------------------------------------------------------------
    % Perform type-2 NUFFT (nonuniform, image space <= uniform, k-space)
    % Siemens: k-space <=> image space
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s:(PCG=%d)(ell=%2d/%2d) Performing type-2 NUFFT... ', datetime, cg_iter, ell, L);
    FhRhUhd = complex(zeros(N, Nc, 'single'));
    cj = type2_scale_factor * finufft2d2(p1, p2, iflag, eps, RhUhd); % N x Nc
    FhRhUhd((support_mask > 0),:) = cj;
    FhRhUhd = reshape(FhRhUhd, [Nkx Nky Nkz Nc]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Calculate (I_{Nc} kron V_{ell}) * (I_{Nc} kron F^H * R_{Omega}^H) * (I_{Nc} kron U_{ell}^H) * d
    %----------------------------------------------------------------------
    VFhRhUhd = VFhRhUhd + bsxfun(@times, V(:,:,:,ell), FhRhUhd);
end

%% Calculate S^H * sum_{ell=1}^{L} (I_{Nc} kron V_{ell}) * (I_{Nc} kron F^H * R_{Omega}^H) * (I_{Nc} kron U_{ell}^H) * d
outp = sum(bsxfun(@times, conj(sens), VFhRhUhd), 4); % Nkx x Nky x Nkz

%% Calculate P^H * S^H * sum_{ell=1}^{L} (I_{Nc} kron V_{ell}) * (I_{Nc} kron F^H * R_{Omega}^H) * (I_{Nc} kron U_{ell}^H) * d
outp = conj(initial_phase) .* outp; % Nkx x Nky x Nkz

%% Vectorize the output
outp = reshape(outp, [N 1]); % Nkx x Nky x Nkz => N x 1

%% Tikhonov regularization
outp = outp + lambda * inp;

end