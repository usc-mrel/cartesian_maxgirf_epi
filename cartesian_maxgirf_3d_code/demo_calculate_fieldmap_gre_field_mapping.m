% demo_calculate_fieldmap_gre_field_mapping.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/20/2025, Last modified: 01/20/2025

%% Clean slate
close all; clearvars -except json_number nr_json_files json_files json_file fieldmap_json_file grad_file_path; clc;

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
%--------------------------------------------------------------------------
% x (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('x_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
x = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% y (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('y_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
y = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% z (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('z_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
z = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Remove readout oversampling
idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(Nkx/2) + 1;
idx2_range = (-floor(Ny/2):ceil(Ny/2)-1).' + floor(Nky/2) + 1;

x_fov = x(idx1_range, idx2_range, :);
y_fov = y(idx1_range, idx2_range, :);
z_fov = z(idx1_range, idx2_range, :);

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
lambda              = json_fieldmap.recon_parameters.lambda;              % l2 regularization parameter
tol                 = json_fieldmap.recon_parameters.tol;                 % PCG tolerance
maxiter             = json_fieldmap.recon_parameters.maxiter;             % PCG maximum iteration 
slice_type          = json_fieldmap.recon_parameters.slice_type;          % type of an excitation slice: "curved" vs "flat"
gnl_correction_flag = json_fieldmap.recon_parameters.gnl_correction_flag; % 1=yes, 0=no

%--------------------------------------------------------------------------
% Number of slices
%--------------------------------------------------------------------------
if isfield(json_fieldmap, 'nr_slices')
    nr_slices = json_fieldmap.nr_slices;
else
    nr_slices = 1;
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
img_fieldmap = zeros([Nkx_fieldmap Nky_fieldmap nr_slices 1 1 2], 'single');  % N1 x N2 x N3 x 1 x 1 x Ne

x_fieldmap = zeros(Nkx_fieldmap, Nky_fieldmap, nr_slices, 'single');
y_fieldmap = zeros(Nkx_fieldmap, Nky_fieldmap, nr_slices, 'single');
z_fieldmap = zeros(Nkx_fieldmap, Nky_fieldmap, nr_slices, 'single');

slice_acq_order = zeros(nr_slices, 1, 'double');

for slice_number = 1:nr_slices
    %----------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
    %----------------------------------------------------------------------
    if nr_slices > 1 % multi-slice
        if mod(nr_slices,2) == 0 % even
            offset1 = 0;
            offset2 = 1;
        else % odd
            offset1 = 1;
            offset2 = 0;
        end
        if slice_number <= ceil(nr_slices / 2)
            actual_slice_number = 2 * slice_number - offset1;
        else
            actual_slice_number = 2 * (slice_number - ceil(nr_slices / 2)) - offset2;
        end
    else
        actual_slice_number = slice_number;
    end
    slice_acq_order(slice_number) = slice_number;

    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('img_type1_slc%d_eco1_gnl%d_%s_i%d_l%4.2f', slice_number, gnl_correction_flag, slice_type, maxiter, lambda);
    cfl_file = fullfile(output_path_fieldmap, img_filename);
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    img_fieldmap(:,:,actual_slice_number,1,1,1) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('img_type1_slc%d_eco2_gnl%d_%s_i%d_l%4.2f', slice_number, gnl_correction_flag, slice_type, maxiter, lambda);
    cfl_file = fullfile(output_path_fieldmap, img_filename);
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    img_fieldmap(:,:,actual_slice_number,1,1,2) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % x (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path_fieldmap, sprintf('x_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    x_fieldmap(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path_fieldmap, sprintf('y_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    y_fieldmap(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path_fieldmap, sprintf('z_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    z_fieldmap(:,:,actual_slice_number) = readcfl(cfl_file);
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

%% Zeropad in k-space (fieldmap)
idx1_range = (-floor(Nkx_fieldmap/2):ceil(Nkx_fieldmap/2)-1).' + floor(2 * Nkx_fieldmap / 2) + 1;
idx2_range = (-floor(Nky_fieldmap/2):ceil(Nky_fieldmap/2)-1).' + floor(2 * Nky_fieldmap / 2) + 1;
idx3_range = (-floor(nr_slices/2):ceil(nr_slices/2)-1).' + floor(2 * nr_slices / 2) + 1;

ksp_fieldmap_zpad = complex(zeros([2 * Nkx_fieldmap 2 * Nky_fieldmap 2 * nr_slices 1 1 2], 'single'));
for echo_number = 1:2
    %----------------------------------------------------------------------
    % Siemens: k-space <=> image space
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s: Applying inverse FFT to move from image space to k-space... ', datetime);
    ksp_fieldmap = img_fieldmap(:,:,:,1,1,echo_number);
    for dim = 1:3
        ksp_fieldmap = sqrt(size(ksp_fieldmap,dim)) * fftshift(ifft(ifftshift(ksp_fieldmap, dim), [], dim), dim);
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    ksp_fieldmap_zpad(idx1_range, idx2_range, idx3_range, :, :, echo_number) = ksp_fieldmap;
end


%% Perform sinc interpolation in image space
%--------------------------------------------------------------------------
% Siemens: k-space <=> image space
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Applying forward FFT to move from k-space to image space... ', datetime);
img_fieldmap_interp = ksp_fieldmap_zpad;
for dim = 1:3
    img_fieldmap_interp = 1 / sqrt(size(img_fieldmap_interp,dim)) * fftshift(fft(ifftshift(img_fieldmap_interp, dim), [], dim), dim);
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Remove readout oversampling
idx1_range = (-floor(2*Nx_fieldmap/2):ceil(2*Nx_fieldmap/2)-1).' + floor(2 * Nkx_fieldmap / 2) + 1;
idx2_range = (-floor(2*Ny_fieldmap/2):ceil(2*Ny_fieldmap/2)-1).' + floor(2 * Nky_fieldmap / 2) + 1;
img_fieldmap_interp = img_fieldmap_interp(idx1_range,idx2_range,:,:,:,:);

%% Perform phase unwrapping
tstart = tic; fprintf('%s: Performing phase unwrapping... ', datetime);
img_phase = angle(img_fieldmap_interp); % [rad]
img_phase = unwrap(img_phase, [], 6); % Nx x Ny x nr_slices x 1 x 1 x 2
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
yy = reshape(permute(img_phase, [6 1 2 3 4 5]), [2 Nx * Ny * Nz]); % [rad]
A = cat(2, TE, ones(2,1));
tstart = tic; fprintf('%s: Performing linear least-squares fitting... ', datetime);
x_ls = inv(A.' * A) * (A.' * yy); % 2 x Nx * Ny * Nz
fieldmap_raw = reshape(x_ls(1,:) / (2 * pi), [Nx Ny Nz]); % [Hz]
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate the regularization term for the B0 field inhomogeneity
weights = zeros(Nx, Ny, Nz, 'single');
for idx3 = 1:Nz
    for idx2 = 1:Ny
        for idx1 = 1:Nx
            %--------------------------------------------------------------
            % Calculate the k-space weight for B0 field inhomogeneity
            %--------------------------------------------------------------
            kx = (-floor(Nx/2) + idx1 - 1) / Nx;
            ky = (-floor(Ny/2) + idx2 - 1) / Ny;
            kz = (-floor(Nz/2) + idx3 - 1) / Nz;
            weights(idx1,idx2,idx3) = 1 / (1 + w_fB0 * (kx^2 + ky^2 + kz^2))^h_fB0;
        end
    end
end

%% Calculate k-space data
%--------------------------------------------------------------------------
% Siemens: k-space <=> image space
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Applying inverse FFT to move from image space to k-space... ', datetime);
ksp = fieldmap_raw;
for dim = 1:3
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
for dim = 1:3
    fieldmap_smooth = 1 / sqrt(size(fieldmap_smooth,dim)) * fftshift(fft(ifftshift(fieldmap_smooth, dim), [], dim), dim);
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
fieldmap_smooth = real(fieldmap_smooth);

%% Calculate a threshold mask
img_abs = abs(img_fieldmap_interp(:,:,:,1,1,1));

mask = (img_abs > mean(img_abs(:)) * 0.20);
mask = bwareaopen(mask, 60); % Keep only blobs with an area of 60 pixels or more.
se = strel('disk', 20, 0);
mask = imclose(mask, se);
se = strel('disk',5);
mask = imdilate(mask,se);

%% Apply a mask
fieldmap_smooth = fieldmap_smooth .* mask;

%% Write a .cfl file
%--------------------------------------------------------------------------
% fieldmap_raw (Nx x Ny x Nz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('fieldmap_raw_gnl%d_%s', gnl_correction_flag, slice_type));
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, fieldmap_raw);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% fieldmap_smooth (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('fieldmap_smooth_gnl%d_%s', gnl_correction_flag, slice_type));
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, fieldmap_smooth);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
