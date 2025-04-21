% demo_step4_cartesian_maxgirf_3d_estimate_csm.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/18/2025, Last modified: 01/18/2025

%% Clean slate
close all; clearvars -except json_number nr_json_files json_files json_file grad_file_path; clc;

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
    ismrmrd_data_file = strrep(json.ismrmrd_data_file, '/', '\');
    output_path       = strrep(json.output_path, '/', '\');
else
    ismrmrd_data_file = json.ismrmrd_data_file;
    output_path       = json.output_path;
end

%--------------------------------------------------------------------------
% Define the BART directory
%--------------------------------------------------------------------------
bart_path = json.bart_path;

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
Lmax          = json.recon_parameters.Lmax;           % maximum rank of the SVD approximation of a higher-order encoding matrix
L             = json.recon_parameters.L;              % rank of the SVD approximation of a higher-order encoding matrix
lambda        = json.recon_parameters.lambda;         % l2 regularization parameter
tol           = json.recon_parameters.tol;            % PCG tolerance
maxiter       = json.recon_parameters.maxiter;        % PCG maximum iteration 
slice_type    = json.recon_parameters.slice_type;     % type of an excitation slice: "curved" vs "flat"
cal_size      = json.recon_parameters.cal_size.';     % size of calibration region
phc_flag      = json.recon_parameters.phc_flag;       % 1=yes, 0=no
gridding_flag = json.recon_parameters.gridding_flag;  % 1=yes, 0=no
cfc_flag      = json.recon_parameters.cfc_flag;       % 1=yes, 0=no
sfc_flag      = json.recon_parameters.sfc_flag;       % 1=yes, 0=no
gnc_flag      = json.recon_parameters.gnc_flag;       % 1=yes, 0=no

if ~cfc_flag
    L = 1;
end

%% Make an output path
mkdir(output_path);

%% Set up BART commands
%--------------------------------------------------------------------------
% Define a BART command
%--------------------------------------------------------------------------
if ispc
    command_prefix = 'wsl';
else
    command_prefix = '';
end
bart_command = sprintf('%s %s/bart', command_prefix, bart_path);

%--------------------------------------------------------------------------
% Translate from a Windows path to a WSL path 
%--------------------------------------------------------------------------
if ispc
    bart_output_path = strrep(output_path, '\', '/');
    bart_output_path = sprintf('/mnt/%s/%s/', lower(bart_output_path(1)), bart_output_path(4:end));
else
    bart_output_path = sprintf('%s/', output_path);
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
encoded_fov(1) = header.encoding.encodedSpace.fieldOfView_mm.x * 1e-3; % [m] RO
encoded_fov(2) = header.encoding.encodedSpace.fieldOfView_mm.y * 1e-3; % [m] PE
encoded_fov(3) = header.encoding.encodedSpace.fieldOfView_mm.z * 1e-3; % [m] SL

Nkx = header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky = header.encoding.encodedSpace.matrixSize.y; % number of phase-encoding steps in k-space
Nkz = header.encoding.encodedSpace.matrixSize.z; % number of slice-encoding steps in k-space

encoded_resolution = encoded_fov ./ [Nkx Nky Nkz]; % [m]

%--------------------------------------------------------------------------
% Number of receive channels
%--------------------------------------------------------------------------
Nc = header.acquisitionSystemInformation.receiverChannels;

%% Read a .cfl file
%--------------------------------------------------------------------------
% ksp_cal_cartesian (Nkx x Nky x Nkz x Nc)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('ksp_cal_cartesian_gridding%d_phc%d', gridding_flag, phc_flag));
if ~exist(strcat(cfl_file, '.cfl'), 'file')
    cfl_file = fullfile(output_path, sprintf('ksp_img_cartesian_gridding%d_phc%d', gridding_flag, phc_flag));
end
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
ksp = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% mask (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'mask_cal');
if ~exist(strcat(cfl_file, '.cfl'), 'file')
    cfl_file = fullfile(output_path, 'mask_img');
end
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
mask = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% circle_mask (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('circle_mask_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
circle_mask = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

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

%--------------------------------------------------------------------------
% dx (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('dx_%s', slice_type));
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
dx = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dy (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('dy_%s', slice_type));
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
dy = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dz (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('dz_%s', slice_type));
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
dz = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% u (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('u_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
u = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% v (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('v_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
v = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% w (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('w_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
w = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% du (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('du_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
du = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dv (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('dv_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
dv = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dw (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('dw_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
dw = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% U (NkNs x Lmax)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('U_cal_%s_cfc%d_sfc%d', slice_type, cfc_flag, sfc_flag));
if ~exist(strcat(cfl_file, '.cfl'), 'file')
    cfl_file = fullfile(output_path, sprintf('U_img_%s_cfc%d_sfc%d', slice_type, cfc_flag, sfc_flag));
end
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
U = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% S (Lmax x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('S_cal_%s_cfc%d_sfc%d', slice_type, cfc_flag, sfc_flag));
if ~exist(strcat(cfl_file, '.cfl'), 'file')
    cfl_file = fullfile(output_path, sprintf('S_img_%s_cfc%d_sfc%d', slice_type, cfc_flag, sfc_flag));
end
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
S = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% V (Nkx x Nky x Nkz x Lmax)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('V_cal_%s_cfc%d_sfc%d', slice_type, cfc_flag, sfc_flag));
if ~exist(strcat(cfl_file, '.cfl'), 'file')
    cfl_file = fullfile(output_path, sprintf('V_img_%s_cfc%d_sfc%d', slice_type, cfc_flag, sfc_flag));
end
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
V = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate parameters for Type-1 and Type-2 NUFFTs
if gnc_flag
    p1 = (u + du) / encoded_fov(1) * (2 * pi); % RO [-pi,pi]
    p2 = (v + dv) / encoded_fov(2) * (2 * pi); % PE [-pi,pi]
    p3 = (w + dw) / encoded_fov(3) * (2 * pi); % SL [-pi,pi]
else
    p1 = u / encoded_fov(1) * (2 * pi); % RO [-pi,pi]
    p2 = v / encoded_fov(2) * (2 * pi); % PE [-pi,pi]
    p3 = w / encoded_fov(3) * (2 * pi); % SL [-pi,pi]
end

%% Calculate a support mask
support_mask = zeros(Nkx, Nky, Nkz, 'single');
support_mask((abs(p1) < pi) & (abs(p2) < pi) & (abs(p3) < pi) & (circle_mask > 0)) = 1;

%% Select spatial positions within a support mask
p1 = p1(support_mask > 0);
p2 = p2(support_mask > 0);
p3 = p3(support_mask > 0);

%% Set parameters for Type-1 and Type-2 NUFFTs
eps = 1e-6;
iflag = -1;

%% Calculate "fake" coil sensitivity maps (Nkx x Nky x Nkz)
sens = complex(ones(Nkx, Nky, Nkz, 'single'));

%% Type-1 NUFFT based Cartesian MaxGIRF operators
Ah = @(x) cartesian_maxgirf_3d_adjoint(x, sens, mask, p1, p2, p3, iflag, eps, support_mask, U, V, L);
AhA = @(x) cartesian_maxgirf_3d_normal(x, sens, mask, p1, p2, p3, iflag, eps, support_mask, U, V, L, lambda);

%% Perform Cartesian MaxGIRF reconstruction
cimg = complex(zeros(Nkx, Nky, Nkz, Nc, 'single'));
for c = 1:Nc
    clear cartesian_maxgirf_3d_normal;
    tstart = tic; fprintf('%s:(c=%2d/%2d) Performing Cartesian MaxGIRF reconstruction:\n', datetime, c, Nc);
    [img, flag, relres, iter, resvec] = pcg(@(x) AhA(x), Ah(ksp(:,:,:,c)), tol, maxiter);
    cimg(:,:,:,c) = reshape(img, [Nkx Nky Nkz]);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
end

%% Write a .cfl file
%--------------------------------------------------------------------------
% cimg (Nkx x Nky x Nkz x Nc)
%--------------------------------------------------------------------------
cimg_filename = sprintf('cimg_cal_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_i%d_l%4.2f', slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, maxiter, lambda);
cfl_file = fullfile(output_path, cimg_filename);
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, cimg);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate gridded k-space
%--------------------------------------------------------------------------
% Siemens: k-space <=> image space
% BART:                image space <=> k-space
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Applying forward FFT to move from image space to k-space... ', datetime);
kgrid = cimg;
for dim = 1:3
    kgrid = 1 / sqrt(size(kgrid,dim)) * fftshift(fft(ifftshift(kgrid, dim), [], dim), dim);
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Write a .cfl file
%--------------------------------------------------------------------------
% kgrid (Nkx x Nky x Nkz x Nc)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('kgrid_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d', slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag));
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, kgrid);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Estimate coil sensitivity using ESPIRiT calibration
%--------------------------------------------------------------------------
% ESPIRiT calibration
%--------------------------------------------------------------------------
% Usage: ecalib [-t f] [-c f] [-k d:d:d] [-r d:d:d] [-m d] [-S] [-W] [-I] [-1]
%               [-P] [-v f] [-a] [-d d] <kspace> <sensitivities> [<ev-maps>]
%
% Estimate coil sensitivities using ESPIRiT calibration.
% Optionally outputs the eigenvalue maps.
%
% -t threshold     This determines the size of the null-space
% -c crop_value    Crop the sensitivities if the eigenvalue is smaller than {crop_value}
% -k ksize         kernel size
% -r cal_size      Limits the size of the calibration region
% -m maps          Number of maps to compute
% -S               create maps with smooth transitions (Soft-SENSE)
% -W               soft-weighting of the singular vectors
% -I               intensity correction
% -1               perform only first part of the calibration
% -P               Do not rotate the phase with respect to the first principal component
% -v variance      Variance of noise in data
% -a               Automatically pick thresholds
% -d level         Debug level
% -h               help
%--------------------------------------------------------------------------
kgrid_file   = strcat(bart_output_path, sprintf('kgrid_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d', slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag));
sens_file    = strcat(bart_output_path, sprintf('sens_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d', slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag));
ev_maps_file = strcat(bart_output_path, sprintf('ev_maps_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d', slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag));
command = sprintf('%s ecalib -t 0.001 -c 0 -k6:6:6 -r%d:%d:%d -m 1 -d5 %s %s %s', bart_command, cal_size(1), cal_size(2), cal_size(3), kgrid_file, sens_file, ev_maps_file);
tstart = tic; fprintf('%s:[BART] Estimating coil sensitivities using ESPIRiT calibration:\n%s\n', datetime, command);
[status_ecalib,result_ecalib] = system(command);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Read a .cfl file
%--------------------------------------------------------------------------
% sens (Nkx x Nky x Nkz x Nc)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('sens_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d', slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
sens = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .cfl file
%--------------------------------------------------------------------------
% ev_maps (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('ev_maps_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d', slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
ev_maps = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Display coil images
slice_number = floor(Nkz/2) + 1;
nr_rows = 2;
nr_cols = 7;
cimg_montage = complex(zeros(Nkx * nr_rows, Nky * nr_cols, 'single'));
sens_montage = complex(zeros(Nkx * nr_rows, Nky * nr_cols, 'single'));

count = 0;
for idx1 = 1:nr_rows
    for idx2 = 1:nr_cols
        idx1_range = (1:Nkx).' + (idx1 - 1) * Nkx;
        idx2_range = (1:Nky).' + (idx2 - 1) * Nky;
        count = count + 1;
        cimg_montage(idx1_range,idx2_range) = cimg(:,:,slice_number,count);
        sens_montage(idx1_range,idx2_range) = sens(:,:,slice_number,count);
        if count >= Nc
            break;
        end
    end
end

figure('Color', 'k', 'Position', [1 5 1239 973]);
imagesc(abs(cimg_montage));
axis image off;
colormap(gray(256));
caxis([0 5]);

fig_filename = sprintf('cimg_cal_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_i%d_l%4.2f', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, maxiter, lambda);

title({'Magnitude of coil images (calibration data)', ...
    sprintf('Cartesian MaxGIRF, SLC = %d, %s slice', slice_number, slice_type), ...
    sprintf('Gridding/PHC/CFC/SFC/GNC = %d/%d/%d/%d/%d', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag)}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
export_fig(fullfile(output_path, sprintf('%s_mag', fig_filename)), '-r300', '-tif', '-c[480, 140, 800, 440]'); % [top,right,bottom,left]
close gcf;

figure('Color', 'k', 'Position', [1 5 1239 973]);
imagesc(angle(cimg_montage) * 180 / pi);
axis image off;
caxis([-180 180]);
colormap(hsv(256));
hc = colorbar;
set(hc, 'Color', 'w', 'FontSize', 14, 'Position', [0.9152 0.2857 0.0121 0.4470], 'TickLabelInterpreter', 'latex');
title(hc, '[deg]', 'Color', 'w', 'Interpreter', 'latex');
title({'Phase of coil images (calibration data)', ...
    sprintf('Cartesian MaxGIRF, SLC = %d, %s slice', slice_number, slice_type), ...
    sprintf('Gridding/PHC/CFC/SFC/GNC = %d/%d/%d/%d/%d', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag)}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
export_fig(fullfile(output_path, sprintf('%s_phase', fig_filename)), '-r300', '-tif', '-c[480, 140, 800, 440]'); % [top,right,bottom,left]
close gcf;

fig_filename = sprintf('sens_cal_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_i%d_l%4.2f', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, maxiter, lambda);

figure('Color', 'k', 'Position', [1 5 1239 973]);
imagesc(abs(sens_montage));
axis image off;
colormap(gray(256));
caxis([0 1]);
title({'Magnitude of coil sensitivity maps', ...
    sprintf('ESPIRiT, SLC = %d, %s slice', slice_number, slice_type), ...
    sprintf('Gridding/PHC/CFC/SFC/GNC = %d/%d/%d/%d/%d', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag)}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
export_fig(fullfile(output_path, sprintf('%s_mag', fig_filename)), '-r300', '-tif', '-c[480, 140, 800, 440]'); % [top,right,bottom,left]
close gcf;

figure('Color', 'k', 'Position', [1 5 1239 973]);
imagesc(angle(sens_montage) * 180 / pi);
axis image off;
caxis([-180 180]);
colormap(hsv(256));
hc = colorbar;
set(hc, 'Color', 'w', 'FontSize', 14, 'Position', [0.9152 0.2857 0.0121 0.4470], 'TickLabelInterpreter', 'latex');
title(hc, '[deg]', 'Color', 'w', 'Interpreter', 'latex');
title({'Phase of coil sensitivity maps', ...
    sprintf('ESPIRiT, SLC = %d, %s slice', slice_number, slice_type), ...
    sprintf('Gridding/PHC/CFC/SFC/GNC = %d/%d/%d/%d/%d', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag)}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
export_fig(fullfile(output_path, sprintf('%s_phase', fig_filename)), '-r300', '-tif', '-c[480, 140, 800, 440]');
close gcf;
