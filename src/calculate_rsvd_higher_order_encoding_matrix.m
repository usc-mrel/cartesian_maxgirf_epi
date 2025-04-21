function [U,S,V] = calculate_rsvd_higher_order_encoding_matrix(k, p, L, s, fieldmap, t, static_B0_correction)
% k: Nk x Nl
% p: N x Nl
% L: rank of approximation
% s: oversampling parameter
% B0map: N x 1, static off-resonance map [Hz]
% t: Nk x 1, time vector [sec]

%% Start a stopwatch timer
start_time = tic;

%% Get dimension parameters
[Nk,Nl] = size(k);
[N,Nl] = size(p);

%% Initialization
W = complex(zeros(Nk, L + s, 'like', k));

%% Set the number of chunks
chunk_size = 2000;
nr_chunks = ceil(N / chunk_size);

%% Calculate W
for idx = 1:nr_chunks
    tstart = tic; fprintf('%s: Calculating W (chunk=%d/%d)... ', datetime, idx, nr_chunks);
    %----------------------------------------------------------------------
    % Calculate the range of columns
    %----------------------------------------------------------------------
    if idx < nr_chunks
        col_range = (1:chunk_size).' + chunk_size * (idx - 1);
    else
        col_range = (1:N - chunk_size * (nr_chunks - 1)).' + chunk_size * (nr_chunks - 1);
    end
    N_chunk = length(col_range);

    %----------------------------------------------------------------------
    % Calculate the columns of E
    %----------------------------------------------------------------------
    phi = k * p(col_range,:).'; % (Nk x Nl) * (N_chunk x Nl).' => Nk x N_chunk
    if static_B0_correction
        % (Nk x 1) * (N_chunk x 1).' => Nk x N_chunk
        % [2pi rad/cycle] * [sec] * [Hz] => [rad]
        phi = phi + (2 * pi) * t * fieldmap(col_range).';
    end
    E_chunk = exp(1j * phi); % Nk x N_chunk

    %----------------------------------------------------------------------
    % Calculate the rows of a Gaussian random matrix, Omega
    %----------------------------------------------------------------------
    Omega_chunk = randn(N_chunk, L + s, 'like', k); % N_chunk x (L + s)

    %----------------------------------------------------------------------
    % Calculate W
    %----------------------------------------------------------------------
    W = W + E_chunk * Omega_chunk; % Nk x (L + s)
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Calculate Y (Nk x (L + s))
Y = W;

%% Calculate the QR decomposition of Y
tstart = tic; fprintf('%s: Calculating the QR decomposition of Y... ', datetime);
[Q,R] = qr(Y,0); % Q:Nk x (L + s), R: (L + s) x (L + s)
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Clear variables
clear W Y;

%% Initialization
B = complex(zeros(L + s, N, 'like', k));

%% Calculate B ((L + s) x N)
for idx = 1:nr_chunks
    tstart = tic; fprintf('%s: Calculating B (chunk=%d/%d)... ', datetime, idx, nr_chunks);
    %----------------------------------------------------------------------
    % Calculate the range of columns
    %----------------------------------------------------------------------
    if idx < nr_chunks
        col_range = (1:chunk_size).' + chunk_size * (idx - 1);
    else
        col_range = (1:N - chunk_size * (nr_chunks - 1)).' + chunk_size * (nr_chunks - 1);
    end

    %----------------------------------------------------------------------
    % Calculate the columns of E
    %----------------------------------------------------------------------
    phi = k * p(col_range,:).'; % (Nk x Nl) * (N_chunk x Nl).' => Nk x N_chunk
    if static_B0_correction
        % (Nk x 1) * (N_chunk x 1).' => Nk x N_chunk
        phi = phi + (2 * pi) * t * fieldmap(col_range).';
    end
    E_chunk = exp(1j * phi); % Nk x N_chunk

    %----------------------------------------------------------------------
    % Calculate the columns of B
    %----------------------------------------------------------------------
    B(:,col_range) = Q' * E_chunk;
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Calculate the SVD of B
tstart = tic; fprintf('%s: Calculating the SVD of B... ', datetime);
[U_hat,S,V] = svds(B,L + s); % U_hat:(L + s) x (L + s), S:(L + s) x (L + s), V:N x (L + s)
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate U
tstart = tic; fprintf('%s: Calculating U... ', datetime);
U = Q * U_hat; % (Nk x L + s) * (L + s x L + s) => Nk x (L + s)
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

end