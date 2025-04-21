function outp = cartesian_maxgirf_adjoint_phase_constrained_lowrank(inp, sens, mask, p1, p2, iflag, eps, support_mask, U, V, L, initial_phase, phase_lowres)

%% Start a stopwatch timer
start_time = tic;

%% Get imaging parameters
[Nkx,Nky,Nkz,Nc] = size(sens);

%--------------------------------------------------------------------------
% Calculate the total number of voxels
%--------------------------------------------------------------------------
N = Nkx * Nky * Nkz;

%--------------------------------------------------------------------------
% Calculate the number of repetitions (averages)
%--------------------------------------------------------------------------
Na = size(phase_lowres,5);

%% Calculate scale factors
type2_scale_factor = 1 / sqrt(Nkx * Nky * Nkz);

%% Calculate D (Nk x Nc x Na)
inp = reshape(inp, [N Nc Na]);
D = inp((mask > 0),:,:); % Nk x Nc x Na

%% Calculate E^H * D = sum_{ell=1}^{L} (I_{Nc} kron V_{ell}) * (I_{Nc} kron F^H * R_{Omega}^H) * (I_{Nc} kron U_{ell}^H) * D
EhD = complex(zeros(Nkx, Nky, Nkz, Nc, Na, 'single'));

for ell = 1:L
    %----------------------------------------------------------------------
    % Calculate (I_{Nc} kron U_{ell}^H) * D
    %----------------------------------------------------------------------
    UhD = bsxfun(@times, conj(U(:,ell)), D); % Nk x Nc x Na

    %----------------------------------------------------------------------
    % Calculate (I_{Nc} kron F^H * R_{Omega}^H) * (I_{Nc} kron U_{ell}^H) * D
    %----------------------------------------------------------------------
    RhUhD = complex(zeros(N, Nc, Na, 'single'));
    RhUhD((mask > 0),:,:) = UhD; % Nk x Nc x Na => N x Nc x Na
    RhUhD = reshape(RhUhD, [Nkx Nky Nkz Nc Na]);

    %----------------------------------------------------------------------
    % Perform type-2 NUFFT (nonuniform, image space <= uniform, k-space)
    % Siemens: k-space <=> image space
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s:(ell=%2d/%2d) Performing type-2 NUFFT... ', datetime, ell, L);
    FhRhUhD = complex(zeros(N, Nc * Na, 'single'));
    cj = type2_scale_factor * finufft2d2(p1, p2, iflag, eps, reshape(RhUhD, [Nkx Nky Nc * Na])); % N x Nc * Na
    FhRhUhD((support_mask > 0),:) = cj;
    FhRhUhD = reshape(FhRhUhD, [Nkx Nky Nkz Nc Na]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Calculate (I_{Nc} kron V_{ell}) * (I_{Nc} kron F^H * R_{Omega}^H) * (I_{Nc} kron U_{ell}^H) * D
    %----------------------------------------------------------------------
    EhD = EhD + bsxfun(@times, V(:,:,:,ell), FhRhUhD);
end

%% Calculate S^H * E^H * D
ShEhD = sum(bsxfun(@times, conj(sens), EhD), 4); % Nkx x Nky x Nkz x 1 x Na

%% Calculate conj(B) \odot S^H * E^H * D
BShEhD = bsxfun(@times, conj(phase_lowres), ShEhD); % Nkx x Nky x Nkz x 1 x Na

%% Calculate P^H * (conj(B) \odot S^H * E^H * D)
PhBShEhD = bsxfun(@times, conj(initial_phase), BShEhD); % Nkx x Nky x Nkz x 1 x Na

%% Calculate P^H * (conj(B) \odot S^H * E^H * D) * Phi^H
outp = sum(PhBShEhD,5); % Nkx x Nky x Nkz

%% Vectorize the output
outp = reshape(outp, [N 1]); % Nk x Nky x Nkz => N x 1

end