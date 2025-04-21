% demo_calculate_fieldmap_gre_field_mapping.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/20/2025, Last modified: 01/20/2025

%% Clean slate
close all; clearvars -except json_number nr_json_files json_files json_file grad_file_path fieldmap_json_file; clc;

%% Start a stopwatch timer
start_time = tic;

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file);
fid = fopen(json_file); 
json_txt = fread(fid, [1 inf], 'char=>char'); 
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
if ispc
    output_path       = strrep(json.output_path, '/', '\');
    ismrmrd_data_file = strrep(json.ismrmrd_data_file, '/', '\');
else
    output_path       = json.output_path;
    ismrmrd_data_file = json.ismrmrd_data_file;
end

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
slice_type = json.recon_parameters.slice_type; % type of an excitation slice: "curved" vs "flat"
sfc_flag   = json.recon_parameters.sfc_flag;   % 1=yes, 0=no

%--------------------------------------------------------------------------
% Number of slices
%--------------------------------------------------------------------------
if isfield(json, 'nr_slices')
    nr_slices = json.nr_slices;
else
    nr_slices = 1;
end

%% Check if static field correction is on
if sfc_flag == 0
    return;
end

%% Read an ISMRMRD file (k-space data)
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... ', datetime, ismrmrd_data_file);
if exist(ismrmrd_data_file, 'file')
    dset = ismrmrd.Dataset(ismrmrd_data_file, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_file);
end

%% Get imaging parameters from an XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% Encoding Space (Nkx, Nky, Nkz)
%--------------------------------------------------------------------------
Nkx = header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky = header.encoding.encodedSpace.matrixSize.y; % number of phase-encoding steps in k-space
Nkz = header.encoding.encodedSpace.matrixSize.z; % number of slice-encoding steps in k-space

%--------------------------------------------------------------------------
% Recon Space (Nx, Ny, Nz)
%--------------------------------------------------------------------------
Nx = header.encoding.reconSpace.matrixSize.x; % number of samples in image space (RO)
Ny = header.encoding.reconSpace.matrixSize.y; % number of samples in image space (PE)
Nz = header.encoding.reconSpace.matrixSize.z; % number of samples in image space (SL)

%% Load a .cfi file
x = zeros(Nkx, Nky, nr_slices, 'single');
y = zeros(Nkx, Nky, nr_slices, 'single');
z = zeros(Nkx, Nky, nr_slices, 'single');

for slice_number = 1:nr_slices
    %----------------------------------------------------------------------
    % x (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('x_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    x(:,:,slice_number) = real(readcfl(cfl_file));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('y_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    y(:,:,slice_number) = real(readcfl(cfl_file));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('z_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    z(:,:,slice_number) = real(readcfl(cfl_file));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, fieldmap_json_file);
fid_fieldmap = fopen(fieldmap_json_file); 
json_txt_fieldmap = fread(fid_fieldmap, [1 inf], 'char=>char'); 
fclose(fid_fieldmap);
json_fieldmap = jsondecode(json_txt_fieldmap);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
if ispc
    output_path_fieldmap       = strrep(json_fieldmap.output_path, '/', '\');
    ismrmrd_data_file_fieldmap = strrep(json_fieldmap.ismrmrd_data_file, '/', '\');
else
    output_path_fieldmap       = json_fieldmap.output_path;
    ismrmrd_data_file_fieldmap = json_fieldmap.ismrmrd_data_file;
end

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
lambda     = json_fieldmap.recon_parameters.lambda;     % l2 regularization parameter
tol        = json_fieldmap.recon_parameters.tol;        % PCG tolerance
maxiter    = json_fieldmap.recon_parameters.maxiter;    % PCG maximum iteration 
slice_type = json_fieldmap.recon_parameters.slice_type; % type of an excitation slice: "curved" vs "flat"
gnc_flag   = json_fieldmap.recon_parameters.gnc_flag;   % 1=yes, 0=no

%--------------------------------------------------------------------------
% Number of slices
%--------------------------------------------------------------------------
if isfield(json_fieldmap, 'nr_slices')
    nr_slices_fieldmap = json_fieldmap.nr_slices;
else
    nr_slices_fieldmap = 1;
end

%--------------------------------------------------------------------------
% Smoothing parameters
%--------------------------------------------------------------------------
if isfield(json_fieldmap, 'smoothing_parameters')
    w_fB0 = json_fieldmap.smoothing_parameters.w_fB0; % 32 in BART?
    h_fB0 = json_fieldmap.smoothing_parameters.h_fB0; % Sobolev index for B0 field inhomogeneity
else
    w_fB0 = 32;
    h_fB0 = 2;
end

h_fB0 = 4;

%% Read an ISMRMRD file (k-space data)
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... ', datetime, ismrmrd_data_file_fieldmap);
if exist(ismrmrd_data_file_fieldmap, 'file')
    dset = ismrmrd.Dataset(ismrmrd_data_file_fieldmap, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_file_fieldmap);
end

%% Get imaging parameters from an XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% Encoding Space (Nkx, Nky, Nkz)
%--------------------------------------------------------------------------
Nkx_fieldmap = header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky_fieldmap = header.encoding.encodedSpace.matrixSize.y; % number of phase-encoding steps in k-space
Nkz_fieldmap = header.encoding.encodedSpace.matrixSize.z; % number of slice-encoding steps in k-space

%--------------------------------------------------------------------------
% Recon Space (Nx, Ny, Nz)
%--------------------------------------------------------------------------
Nx_fieldmap = header.encoding.reconSpace.matrixSize.x; % number of samples in image space (RO)
Ny_fieldmap = header.encoding.reconSpace.matrixSize.y; % number of samples in image space (PE)
Nz_fieldmap = header.encoding.reconSpace.matrixSize.z; % number of samples in image space (SL)

%% Load a .cfl file
%--------------------------------------------------------------------------
% TE (N1 x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path_fieldmap, 'TE');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
TE = real(readcfl(cfl_file)).';
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Load a .cfl file
img_fieldmap = zeros([Nkx_fieldmap Nky_fieldmap nr_slices_fieldmap 1 1 2], 'single');  % N1 x N2 x N3 x 1 x 1 x Ne

x_fieldmap = zeros(Nkx_fieldmap, Nky_fieldmap, nr_slices_fieldmap, 'single');
y_fieldmap = zeros(Nkx_fieldmap, Nky_fieldmap, nr_slices_fieldmap, 'single');
z_fieldmap = zeros(Nkx_fieldmap, Nky_fieldmap, nr_slices_fieldmap, 'single');

for slice_number = 1:nr_slices_fieldmap
    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('img_type1_slc%d_eco1_gnc%d_%s_i%d_l%4.2f', slice_number, gnc_flag, slice_type, maxiter, lambda);
    cfl_file = fullfile(output_path_fieldmap, img_filename);
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    img_fieldmap(:,:,slice_number,1,1,1) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('img_type1_slc%d_eco2_gnc%d_%s_i%d_l%4.2f', slice_number, gnc_flag, slice_type, maxiter, lambda);
    cfl_file = fullfile(output_path_fieldmap, img_filename);
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    img_fieldmap(:,:,slice_number,1,1,2) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % x (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path_fieldmap, sprintf('x_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    x_fieldmap(:,:,slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path_fieldmap, sprintf('y_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    y_fieldmap(:,:,slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path_fieldmap, sprintf('z_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    z_fieldmap(:,:,slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%--------------------------------------------------------------------------
% read_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'read_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
read_sign_fieldmap = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% phase_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'phase_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
phase_sign_fieldmap = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform phase unwrapping
tstart = tic; fprintf('%s: Performing phase unwrapping... ', datetime);
img_phase = angle(img_fieldmap); % [rad]
img_phase = unwrap(img_phase, [], 6); % Nkx_fieldmap x Nky_fieldmap x nr_slices x 1 x 1 x 2
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Normalize the phase
tstart = tic; fprintf('%s: Normalize the unwrapped phase of multi-echo images... ', datetime);
img_phase = img_phase - repmat(img_phase(:,:,:,:,:,1), [1 1 1 1 1 2]);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform a linear least-squares fit of the phase-time curve
%--------------------------------------------------------------------------
% y in [rad], a in [rad/sec], TE in [sec]
% y1 = a * TE1 + b    [y1]   [TE1 1] [a]
% y2 = a * TE2 + b => [y2] = [TE2 1] [b] => y = A * x
% yn = a * TEn + b    [yn]   [TEn 1]
%
% x_ls = inv(A.' * A) * A.' * y
%--------------------------------------------------------------------------
yy = reshape(permute(img_phase, [6 1 2 3 4 5]), [2 Nkx_fieldmap * Nky_fieldmap * nr_slices_fieldmap]); % [rad]
A = cat(2, TE, ones(2,1));
tstart = tic; fprintf('%s: Performing linear least-squares fitting... ', datetime);
x_ls = inv(A.' * A) * (A.' * yy); % 2 x Nkx_fieldmap * Nky_fieldmap * nr_slices_fieldmap
fieldmap_raw = reshape(x_ls(1,:) / (2 * pi), [Nkx_fieldmap Nky_fieldmap nr_slices_fieldmap]); % [Hz]
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform linear interpolation on a fieldmap
tstart = tic; fprintf('%s: Perform linear interpolation on a fieldmap... ', datetime);
F = scatteredInterpolant(double(x_fieldmap(:)), double(y_fieldmap(:)), double(z_fieldmap(:)), double(fieldmap_raw(:)), 'linear', 'nearest');
fieldmap_interp = reshape(F(double(x(:)), double(y(:)), double(z(:))), [Nkx Nky nr_slices]);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate a threshold mask
img_abs = sum(abs(img_fieldmap(:,:,:,1,1,1)),6);
mask = (img_abs > mean(img_abs(:)) * 0.60);
mask = bwareaopen(mask, 60); % Keep only blobs with an area of 60 pixels or more.

%% Perform linear interpolation on a mask
tstart = tic; fprintf('%s: Perform linear interpolation on a mask... ', datetime);
F = scatteredInterpolant(double(x_fieldmap(:)), double(y_fieldmap(:)), double(z_fieldmap(:)), double(mask(:)), 'linear', 'nearest');
mask_interp = reshape(F(double(x(:)), double(y(:)), double(z(:))), [Nkx Nky nr_slices]);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Apply a threshold mask
fieldmap_interp_masked = fieldmap_interp .* mask_interp;

%% Calculate the regularization term for the B0 field inhomogeneity
weights = zeros(Nkx, Nky, 'single');
for idx2 = 1:Nky
    for idx1 = 1:Nkx
        %------------------------------------------------------------------
        % Calculate the k-space weight for B0 field inhomogeneity
        %------------------------------------------------------------------
        kx = (-floor(Nkx/2) + idx1 - 1) / Nkx;
        ky = (-floor(Nky/2) + idx2 - 1) / Nky;
        weights(idx1,idx2) = 1 / (1 + w_fB0 * (kx^2 + ky^2))^h_fB0;
    end
end

%% Calculate k-space data
%--------------------------------------------------------------------------
% Siemens: k-space <=> image space
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Applying inverse FFT to move from image space to k-space... ', datetime);
ksp = fieldmap_interp_masked;
for dim = 1:2
    ksp = sqrt(size(ksp,dim)) * fftshift(ifft(ifftshift(ksp, dim), [], dim), dim);
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform smoothing via k-space filtering
ksp = bsxfun(@times, weights, ksp);

%% Calculate a smoothed fieldmap [Hz]
%--------------------------------------------------------------------------
% Siemens: k-space <=> image space
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Applying forward FFT to move from k-space to image space... ', datetime);
fieldmap_smooth = ksp;
for dim = 1:2
    fieldmap_smooth = 1 / sqrt(size(fieldmap_smooth,dim)) * fftshift(fft(ifftshift(fieldmap_smooth, dim), [], dim), dim);
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
fieldmap_smooth = real(fieldmap_smooth);

%% Display a raw fieldmap
figure('Color', 'w');
montage(fieldmap_raw, 'DisplayRange', []);
colormap(colorcet('D1'));
colorbar;
clim([-15 15]);
export_fig(fullfile(output_path, 'fieldmap_montage_raw'), '-r300', '-tif');

%% Display a interpolated fieldmap
figure('Color', 'w');
montage(fieldmap_interp, 'DisplayRange', []);
colormap(colorcet('D1'));
colorbar;
clim([-15 15]);
export_fig(fullfile(output_path, 'fieldmap_montage_interpolated'), '-r300', '-tif');

%% Display a masked, interpolated fieldmap
figure('Color', 'w');
montage(fieldmap_interp_masked, 'DisplayRange', []);
colormap(colorcet('D1'));
colorbar;
clim([-15 15]);
export_fig(fullfile(output_path, 'fieldmap_montage_masked'), '-r300', '-tif');

%% Display a filtered, masked, interpolated fieldmap
figure('Color', 'w');
montage(fieldmap_smooth, 'DisplayRange', []);
colormap(colorcet('D1'));
colorbar;
clim([-15 15]);
export_fig(fullfile(output_path, 'fieldmap_montage_masked_filtered'), '-r300', '-tif');

%% Write a .cfl file
for slice_number = 1:nr_slices
    %----------------------------------------------------------------------
    % fieldmap_smooth (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('fieldmap_smooth_slc%d_gnc%d_%s', slice_number, gnc_flag, slice_type);
    cfl_file = fullfile(output_path, img_filename);
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, fieldmap_smooth(:,:,slice_number));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end
