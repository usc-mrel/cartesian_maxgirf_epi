function outp = cartesian_maxgirf_3d_adjoint(inp, sens, mask, p1, p2, p3, iflag, eps, support_mask, U, V, L)

%% Start a stopwatch timer
start_time = tic;

%% Get imaging parameters
[Nkx,Nky,Nkz,Nc] = size(sens);

%--------------------------------------------------------------------------
% Calculate the total number of voxels
%--------------------------------------------------------------------------
N = Nkx * Nky * Nkz;

%% Calculate scale factors
type2_scale_factor = 1 / sqrt(Nkx * Nky * Nkz);

%% Calculate d (NkNs x Nc)
inp = reshape(inp, [N Nc]);
d = inp((mask > 0), :); % NkNs x Nc

%% Calculate sum_{ell=1}^{L} (I_{Nc} kron V_{ell}) * (I_{Nc} kron F^H * R_{Omega}^H) * (I_{Nc} kron U_{ell}^H) * d
VFhRhUhd = complex(zeros(Nkx, Nky, Nkz, Nc, 'single'));

for ell = 1:L
    %----------------------------------------------------------------------
    % Calculate (I_{Nc} kron U_{ell}^H) * d
    %----------------------------------------------------------------------
    Uhd = bsxfun(@times, conj(U(:,ell)), d); % NkNs x Nc

    %----------------------------------------------------------------------
    % Calculate (I_{Nc} kron F^H * R_{Omega}^H) * (I_{Nc} kron U_{ell}^H) * d
    %----------------------------------------------------------------------
    RhUhd = complex(zeros(N, Nc, 'single'));
    RhUhd((mask > 0),:) = Uhd; % NkNs x Nc => N x Nc
    RhUhd = reshape(RhUhd, [Nkx Nky Nkz Nc]);

    %----------------------------------------------------------------------
    % Perform type-2 NUFFT (nonuniform, image space <= uniform, k-space)
    % Siemens: k-space <=> image space
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s:(ell=%2d/%2d) Performing type-2 NUFFT... ', datetime, ell, L);
    FhRhUhd = complex(zeros(N, Nc, 'single'));
    cj = type2_scale_factor * finufft3d2(p1, p2, p3, iflag, eps, RhUhd); % N x Nc
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

%% Vectorize the output
outp = reshape(outp, [N 1]); % NkNs x Nky x Nkz => N x 1

end