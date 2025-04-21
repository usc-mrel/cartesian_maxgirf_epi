% demo_step5_cartesian_maxgirf_3d_recon.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/18/2025, Last modified: 01/18/2025

%% Clean slate
close all; clearvars -except json_number nr_json_files json_files json_file grad_file_path; clc;

%% Set a flag to save a figure
save_figure = 1;

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
    siemens_twix_file  = strrep(json.siemens_twix_file, '/', '\');
    ismrmrd_data_file  = strrep(json.ismrmrd_data_file, '/', '\');
    ismrmrd_noise_file = strrep(json.ismrmrd_noise_file, '/', '\');
    output_path        = strrep(json.output_path, '/', '\');
else
    siemens_twix_file  = json.siemens_twix_file;
    ismrmrd_data_file  = json.ismrmrd_data_file;
    ismrmrd_noise_file = json.ismrmrd_noise_file;
    output_path        = json.output_path;
end

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

%--------------------------------------------------------------------------
% main_orientation (SAGITTAL/CORONAL/TRANSVERSAL = 0/1/2)
%--------------------------------------------------------------------------
if isfield(json, 'main_orientation')
    main_orientation = json.main_orientation;
else
    main_orientation = 2;
end

%% Make an output path
mkdir(output_path);

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

%% Read a .cfl file
%--------------------------------------------------------------------------
% ksp_cartesian (Nkx x Nky x Nkz x Nc)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('ksp_img_cartesian_gridding%d_phc%d', gridding_flag, phc_flag));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
ksp = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% mask (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'mask_img');
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
% sens (Nkx x Nky x Nkz x Nc)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('sens_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d', slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
sens = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% ev_maps (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('ev_maps_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d', slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
ev_maps = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% x_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('x_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
x_shift = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% y_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('y_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
y_shift = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% z_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('z_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
z_shift = readcfl(cfl_file);
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
% read_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'read_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
read_sign = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% phase_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'phase_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
phase_sign = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% U (NkNs x Lmax)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('U_img_%s_cfc%d_sfc%d', slice_type, cfc_flag, sfc_flag));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
U = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% S (Lmax x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('S_img_%s_cfc%d_sfc%d', slice_type, cfc_flag, sfc_flag));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
S = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% V (Nkx x Nky x Nkz x Lmax)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('V_img_%s_cfc%d_sfc%d', slice_type, cfc_flag, sfc_flag));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
V = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Apply a threshold mask on CSMs
ev_mask = (ev_maps > 0.94);
ev_mask2 = bwareaopen(ev_mask, 60); % Keep only blobs with an area of 60 pixels or more.
se = strel('disk',5);
ev_mask_dilated = imdilate(ev_mask2,se);
%sens = bsxfun(@times, sens, ev_mask_dilated);

%% Calculate spatial positions for Type-3 NUFFT
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

%% Type-1 NUFFT based Cartesian MaxGIRF operators
Ah = @(x) cartesian_maxgirf_3d_adjoint(x, sens, mask, p1, p2, p3, iflag, eps, support_mask, U, V, L);
AhA = @(x) cartesian_maxgirf_3d_normal(x, sens, mask, p1, p2, p3, iflag, eps, support_mask, U, V, L, lambda);

%% Perform Cartesian MaxGIRF reconstruction
clear cartesian_maxgirf_3d_normal;
tstart = tic; fprintf('%s: Performing Cartesian MaxGIRF reconstruction:\n', datetime);
[img, flag, relres, iter, resvec] = pcg(@(x) AhA(x), Ah(ksp), tol, maxiter);
img = reshape(img, [Nkx Nky Nkz]);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Write a .cfl file
%--------------------------------------------------------------------------
% img (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
img_filename = sprintf('img_maxgirf_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_%s_i%d_l%4.2f', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, slice_type, maxiter, lambda);
cfl_file = fullfile(output_path, img_filename);
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, img);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% support_mask (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
mask_filename = sprintf('support_mask_%s', slice_type);
cfl_file = fullfile(output_path, mask_filename);
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, support_mask);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Display the magnitude of an image
FontSize = 14;

c1 = floor(Nkx/2) + 1;
c2 = floor(Nky/2) + 1;
c3 = floor(Nkz/2) + 1;

slice_number = c3-30;

xmax = max([abs(max(x_shift(:))) abs(min(x_shift(:)))]);
ymax = max([abs(max(y_shift(:))) abs(min(y_shift(:)))]);
zmax = max([abs(max(z_shift(:))) abs(min(z_shift(:)))]);

xlimits = [-xmax xmax];
ylimits = [-ymax ymax];
zlimits = [-zmax zmax];

if main_orientation == 2 % TRANSVERSAL = 2
    Position = [1020 44 525 934];
    slice_direction = 'z';
    slice_offset = z_shift(c1,c2,slice_number) * 1e3;
    ax = 0;
    el = 90;
elseif main_orientation == 0 % SAGITTAL = 0
    Position = [1038 25 778 953];
    slice_direction = 'x';
    slice_offset = x_shift(c1,c2,slice_number) * 1e3;
    ax = -90;
    el = 0;
elseif main_orientation == 1 % CORONAL = 1
    Position = [680 28 1029 950];
    slice_direction = 'y';
    slice_offset = y_shift(c1,c2,slice_number) * 1e3;
    ax = 0;
    el = 0;
end

if read_sign == -1
    x_shift = flip(x_shift,1);
    y_shift = flip(y_shift,1);
    z_shift = flip(z_shift,1);
    img = flip(img,1);
end

if phase_sign == -1
    x_shift = flip(x_shift,2);
    y_shift = flip(y_shift,2);
    z_shift = flip(z_shift,2);
    img = flip(img,2);
end

figure('Color', 'w', 'Position', Position);
surf(x_shift(:,:,slice_number) * 1e3, y_shift(:,:,slice_number) * 1e3, z_shift(:,:,slice_number) * 1e3, abs(img(:,:,slice_number)), 'EdgeColor', 'none');
colormap(gray(256));
axis image;
xlim(xlimits * 1e3);
ylim(ylimits * 1e3);
zlim(zlimits * 1e3);
%caxis([0 15]);
title({'Cartesian MaxGIRF', sprintf('SLC = %d, %s slice, %s = %4.1f mm', slice_number, slice_type, slice_direction, slice_offset), ...
    sprintf('Gridding/PHC/CFC/SFC/GNC = %d/%d/%d/%d/%d', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag), ...
    sprintf('max. iterations = %d, $$\\lambda$$ = %4.2f', maxiter, lambda)}, ...
    'FontWeight', 'normal', 'FontSize', FontSize, 'Interpreter', 'latex');
xlabel('x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
ylabel('y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
zlabel('z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
set(gca, 'TickLabelInterpreter', 'latex');
view(ax,el);
set(gca, 'ZDir', 'reverse');
drawnow;
fig_filename = sprintf('img_maxgirf_slc%d_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_%s_i%d_l%4.2f_mag', slice_number, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, slice_type, maxiter, lambda);
export_fig(fullfile(output_path, fig_filename), '-r300', '-tif');
close gcf;

%% Display the phase of an image
figure('Color', 'w', 'Position', Position);
surf(x_shift(:,:,slice_number) * 1e3, y_shift(:,:,slice_number) * 1e3, z_shift(:,:,slice_number) * 1e3, angle(img(:,:,slice_number)) * 180 / pi, 'EdgeColor', 'none');
colormap(hsv(256));
axis image;
xlim(xlimits * 1e3);
ylim(ylimits * 1e3);
zlim(zlimits * 1e3);
caxis([-180 180]);
title({'Cartesian MaxGIRF', sprintf('SLC = %d, %s slice, %s = %4.1f mm', slice_number, slice_type, slice_direction, slice_offset), ...
    sprintf('Gridding/PHC/CFC/SFC/GNC = %d/%d/%d/%d/%d', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag), ...
    sprintf('max. iterations = %d, $$\\lambda$$ = %4.2f', maxiter, lambda)}, ...
    'FontWeight', 'normal', 'FontSize', FontSize, 'Interpreter', 'latex');
xlabel('x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
ylabel('y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
zlabel('z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
set(gca, 'TickLabelInterpreter', 'latex');
view(ax,el);
set(gca, 'ZDir', 'reverse');
drawnow;
fig_filename = sprintf('img_maxgirf_slc%d_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_%s_i%d_l%4.2f_phase', slice_number, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, slice_type, maxiter, lambda);
export_fig(fullfile(output_path, fig_filename), '-r300', '-tif');
close gcf;

%% Display a support mask
figure('Color', 'w', 'Position', Position);
surf(x_shift(:,:,slice_number) * 1e3, y_shift(:,:,slice_number) * 1e3, z_shift(:,:,slice_number) * 1e3, support_mask(:,:,slice_number), 'EdgeColor', 'none');
colormap(gray(256));
axis image;
xlim(xlimits * 1e3);
ylim(ylimits * 1e3);
zlim(zlimits * 1e3);
caxis([0 1]);
title({'Cartesian MaxGIRF (mask)', sprintf('SLC = %d, %s slice, %s = %4.1f mm', slice_number, slice_type, slice_direction, slice_offset), ...
    sprintf('Gridding/PHC/CFC/SFC/GNC = %d/%d/%d/%d/%d', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag), ...
    sprintf('max. iterations = %d, $$\\lambda$$ = %4.2f', maxiter, lambda)}, ...
    'FontWeight', 'normal', 'FontSize', FontSize, 'Interpreter', 'latex');
xlabel('x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
ylabel('y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
zlabel('z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
set(gca, 'TickLabelInterpreter', 'latex', 'Color', 'k');
view(ax,el);
set(gca, 'ZDir', 'reverse');
drawnow;
export_fig(fullfile(output_path, mask_filename), '-r300', '-tif');
close gcf;
