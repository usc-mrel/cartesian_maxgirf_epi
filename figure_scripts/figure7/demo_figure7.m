% demo_figure7.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 02/17/2025, Last modified: 03/11/2025

%% Clean slate
close all; clear all; clc;

%% Start a stopwatch timer
start_time = tic;

%% Set source directories
package_path = 'D:\cartesian_maxgirf_epi_3d';

%% Add source directories to search path
addpath(genpath(package_path));

%% Define the full path of an output directory
json_file1 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\meas_MID00180_FID08983_ep3d_tra_AP_TR153_TE58_etl61_0_8mm_gridding1_phc1_cfc0_sfc0_gnc0.json';
json_file2 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\meas_MID00180_FID08983_ep3d_tra_AP_TR153_TE58_etl61_0_8mm_gridding1_phc1_cfc0_sfc0_gnc1.json';
json_file3 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\meas_MID00180_FID08983_ep3d_tra_AP_TR153_TE58_etl61_0_8mm_gridding1_phc1_cfc1_sfc0_gnc0.json';
json_file4 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\meas_MID00180_FID08983_ep3d_tra_AP_TR153_TE58_etl61_0_8mm_gridding1_phc1_cfc1_sfc0_gnc1.json';

dicom_path1 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\dicom\S62_ep3d_tra_AP_TR153_TE58_etl61_0.8mm'; % ND
dicom_path2 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\dicom\S63_ep3d_tra_AP_TR153_TE58_etl61_0.8mm_S62_DIS3D'; % DIS3D

%E:\projects_lenovo_20250319\cartesian_maxgirf_epi_2d\data

json_file1 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\meas_MID00180_FID08983_ep3d_tra_AP_TR153_TE58_etl61_0_8mm_gridding1_phc1_cfc0_sfc0_gnc0.json';
json_file2 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\meas_MID00180_FID08983_ep3d_tra_AP_TR153_TE58_etl61_0_8mm_gridding1_phc1_cfc0_sfc0_gnc1.json';
json_file3 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\meas_MID00180_FID08983_ep3d_tra_AP_TR153_TE58_etl61_0_8mm_gridding1_phc1_cfc1_sfc0_gnc0.json';
json_file4 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\meas_MID00180_FID08983_ep3d_tra_AP_TR153_TE58_etl61_0_8mm_gridding1_phc1_cfc1_sfc0_gnc1.json';

dicom_path1 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\dicom\S62_ep3d_tra_AP_TR153_TE58_etl61_0.8mm'; % ND
dicom_path2 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250209\dicom\S63_ep3d_tra_AP_TR153_TE58_etl61_0.8mm_S62_DIS3D'; % DIS3D

%% Get directory information
dir_info = dir(fullfile(dicom_path1, '*IMA'));
nr_files = length(dir_info);

%% Get a DICOM header
dicom_info = dicominfo(fullfile(dir_info(1).folder, dir_info(1).name));

%% Parse the DICOM header
%--------------------------------------------------------------------------
% Rows Attribute
% Number of rows in the image
%--------------------------------------------------------------------------
N1 = double(dicom_info.Rows);

%--------------------------------------------------------------------------
% Columns  Attribute
% Number of columns in the image
%--------------------------------------------------------------------------
N2 = double(dicom_info.Columns);

%% Read a .dicom file
img1_dicom = zeros(N1, N2, nr_files, 'single');
x_dicom1 = zeros(N1, N2, nr_files, 'single');
y_dicom1 = zeros(N1, N2, nr_files, 'single');
z_dicom1 = zeros(N1, N2, nr_files, 'single');

for idx = 1:nr_files

    dicom_file = fullfile(dir_info(idx).folder, dir_info(idx).name);
    tstart = tic; fprintf('%s:(i=%3d/%3d) Reading a .dicom file: %s... ', datetime, idx, nr_files, dicom_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Get a DICOM header
    dicom_info = dicominfo(dicom_file);

    %% Read a DICOM image
    if isfield(dicom_info, 'RescaleSlope')
        RescaleSlope = dicom_info.RescaleSlope;
    else
        RescaleSlope = 1;
    end

    if isfield(dicom_info, 'RescaleIntercept')
        RescaleIntercept = dicom_info.RescaleIntercept;
    else
        RescaleIntercept = 0;
    end

    img1_dicom(:,:,idx) = RescaleSlope * double(dicomread(dicom_info)).' + RescaleIntercept; % transpose it!

    %% Parse the DICOM header
    %----------------------------------------------------------------------
    % Patient Position Attribute
    % Patient position descriptor relative to the equipment
    %----------------------------------------------------------------------
    patient_position = dicom_info.PatientPosition;

    %----------------------------------------------------------------------
    % Slice Thickness Attribute
    % Nominal slice thickness, in mm
    %----------------------------------------------------------------------
    slice_thickness = dicom_info.SliceThickness; % [mm]

    %----------------------------------------------------------------------
    % Image Position (Patient) Attribute
    % The x, y, and z coordinates of the upper left hand corner
    % (center of the first voxel transmitted) of the image, in mm
    %----------------------------------------------------------------------
    ipp = dicom_info.ImagePositionPatient; % [mm]

    %----------------------------------------------------------------------
    % Image Orientation (Patient) Attribute
    % The direction cosines of the first row and the first column with respect
    % to the patient
    %----------------------------------------------------------------------
    iop = dicom_info.ImageOrientationPatient;

    %----------------------------------------------------------------------
    % Pixel Spacing Attribute
    % Physical distance in the Patient between the center of each pixel, specified
    % by a numeric pair - adjacent row spacing, adjacent column spacing in mm
    %----------------------------------------------------------------------
    pixel_spacing = dicom_info.PixelSpacing; % [mm]

    %----------------------------------------------------------------------
    % Number of slices
    %----------------------------------------------------------------------
    N3 = 1;

    %----------------------------------------------------------------------
    % Slice number
    %----------------------------------------------------------------------
    instance_number = dicom_info.InstanceNumber;

    %% Calculate the total number of voxels
    N = N1 * N2 * N3;

    %% Calculate a rotation matrix from the PCS to the DCS
    R_pcs2dcs = siemens_calculate_transform_pcs_to_dcs(patient_position);

    %% Calculate a scaling matrix
    scaling_matrix_dicom = [pixel_spacing(1) 0 0; 0 pixel_spacing(2) 0; 0 0 slice_thickness] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Calculate a trannsformation matrix from the RCS to the PCS [r,c,s] <=> [L,P,S]
    R_rcs2pcs = cat(2, iop(1:3), iop(4:6), cross(iop(1:3), iop(4:6)));

    %% Calculate spatial coordinates in the RCS [m]
    [I1,I2,I3] = ndgrid((0:N1-1).', (0:N2-1).', (0:N3-1).');
    r_rcs = scaling_matrix_dicom * cat(1, I1(:).', I2(:).', I3(:).'); % 3 x N

    %% Calculate spatial coordinates in the LPH [m] (R => L, A => P, I => S)
    % The DICOM LPH coordinate system is identical to the Siemens Patient Coordinate System (PCS)
    r_pcs = repmat(ipp * 1e-3, [1 N]) + R_rcs2pcs * r_rcs;

    %% Add a table position
    r_pcs = r_pcs + repmat(-dicom_info.Private_0019_1014 * 1e-3, [1 N]);

    %% Calculate spatial coordinates in the DCS [m]
    r_dcs = R_pcs2dcs * r_pcs;

    %% Save arrays
    x_dicom1(:,:,idx) = reshape(r_dcs(1,:), [N1 N2]); % N x 1 [m]
    y_dicom1(:,:,idx) = reshape(r_dcs(2,:), [N1 N2]); % N x 1 [m]
    z_dicom1(:,:,idx) = reshape(r_dcs(3,:), [N1 N2]); % N x 1 [m]
end

%% Get directory information
dir_info = dir(fullfile(dicom_path2, '*IMA'));
nr_files = length(dir_info);

%% Get a DICOM header
dicom_info = dicominfo(fullfile(dir_info(1).folder, dir_info(1).name));

%% Parse the DICOM header
%--------------------------------------------------------------------------
% Rows Attribute
% Number of rows in the image
%--------------------------------------------------------------------------
N1 = double(dicom_info.Rows);

%--------------------------------------------------------------------------
% Columns  Attribute
% Number of columns in the image
%--------------------------------------------------------------------------
N2 = double(dicom_info.Columns);

%% Read a .dicom file
img2_dicom = zeros(N1, N2, nr_files, 'single');
x_dicom2 = zeros(N1, N2, nr_files, 'single');
y_dicom2 = zeros(N1, N2, nr_files, 'single');
z_dicom2 = zeros(N1, N2, nr_files, 'single');

for idx = 1:nr_files

    dicom_file = fullfile(dir_info(idx).folder, dir_info(idx).name);
    tstart = tic; fprintf('%s:(i=%3d/%3d) Reading a .dicom file: %s... ', datetime, idx, nr_files, dicom_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Get a DICOM header
    dicom_info = dicominfo(dicom_file);

    %% Read a DICOM image
    if isfield(dicom_info, 'RescaleSlope')
        RescaleSlope = dicom_info.RescaleSlope;
    else
        RescaleSlope = 1;
    end

    if isfield(dicom_info, 'RescaleIntercept')
        RescaleIntercept = dicom_info.RescaleIntercept;
    else
        RescaleIntercept = 0;
    end

    img2_dicom(:,:,idx) = RescaleSlope * double(dicomread(dicom_info)).' + RescaleIntercept; % transpose it!

    %% Parse the DICOM header
    %----------------------------------------------------------------------
    % Patient Position Attribute
    % Patient position descriptor relative to the equipment
    %----------------------------------------------------------------------
    patient_position = dicom_info.PatientPosition;

    %----------------------------------------------------------------------
    % Slice Thickness Attribute
    % Nominal slice thickness, in mm
    %----------------------------------------------------------------------
    slice_thickness = dicom_info.SliceThickness; % [mm]

    %----------------------------------------------------------------------
    % Image Position (Patient) Attribute
    % The x, y, and z coordinates of the upper left hand corner
    % (center of the first voxel transmitted) of the image, in mm
    %----------------------------------------------------------------------
    ipp = dicom_info.ImagePositionPatient; % [mm]

    %----------------------------------------------------------------------
    % Image Orientation (Patient) Attribute
    % The direction cosines of the first row and the first column with respect
    % to the patient
    %----------------------------------------------------------------------
    iop = dicom_info.ImageOrientationPatient;

    %----------------------------------------------------------------------
    % Pixel Spacing Attribute
    % Physical distance in the Patient between the center of each pixel, specified
    % by a numeric pair - adjacent row spacing, adjacent column spacing in mm
    %----------------------------------------------------------------------
    pixel_spacing = dicom_info.PixelSpacing; % [mm]

    %----------------------------------------------------------------------
    % Number of slices
    %----------------------------------------------------------------------
    N3 = 1;

    %----------------------------------------------------------------------
    % Slice number
    %----------------------------------------------------------------------
    instance_number = dicom_info.InstanceNumber;

    %% Calculate the total number of voxels
    N = N1 * N2 * N3;

    %% Calculate a rotation matrix from the PCS to the DCS
    R_pcs2dcs = siemens_calculate_transform_pcs_to_dcs(patient_position);

    %% Calculate a scaling matrix
    scaling_matrix_dicom = [pixel_spacing(1) 0 0; 0 pixel_spacing(2) 0; 0 0 slice_thickness] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Calculate a trannsformation matrix from the RCS to the PCS [r,c,s] <=> [L,P,S]
    R_rcs2pcs = cat(2, iop(1:3), iop(4:6), cross(iop(1:3), iop(4:6)));

    %% Calculate spatial coordinates in the RCS [m]
    [I1,I2,I3] = ndgrid((0:N1-1).', (0:N2-1).', (0:N3-1).');
    r_rcs = scaling_matrix_dicom * cat(1, I1(:).', I2(:).', I3(:).'); % 3 x N

    %% Calculate spatial coordinates in the LPH [m] (R => L, A => P, I => S)
    % The DICOM LPH coordinate system is identical to the Siemens Patient Coordinate System (PCS)
    r_pcs = repmat(ipp * 1e-3, [1 N]) + R_rcs2pcs * r_rcs;

    %% Add a table position
    r_pcs = r_pcs + repmat(-dicom_info.Private_0019_1014 * 1e-3, [1 N]);

    %% Calculate spatial coordinates in the DCS [m]
    r_dcs = R_pcs2dcs * r_pcs;

    %% Save arrays
    x_dicom2(:,:,idx) = reshape(r_dcs(1,:), [N1 N2]); % N x 1 [m]
    y_dicom2(:,:,idx) = reshape(r_dcs(2,:), [N1 N2]); % N x 1 [m]
    z_dicom2(:,:,idx) = reshape(r_dcs(3,:), [N1 N2]); % N x 1 [m]
end

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file1);
fid = fopen(json_file1);
json_txt = fread(fid, [1 inf], 'char=>char');
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
output_path1 = strrep(json.output_path, '/', '\');
% Hack
output_path1 = strrep(output_path1, 'D:', 'E:\projects_lenovo_20250319');

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
Lmax          = json.recon_parameters.Lmax;           % maximum rank of the SVD approximation of a higher-order encoding matrix
L             = json.recon_parameters.L;              % rank of the SVD approximation of a higher-order encoding matrix
lambda        = json.recon_parameters.lambda;         % l2 regularization parameter
tol           = json.recon_parameters.tol;            % PCG tolerance
maxiter       = json.recon_parameters.maxiter;        % PCG maximum iteration 
slice_type    = json.recon_parameters.slice_type;     % type of an excitation slice: "curved" vs "flat"
phc_flag      = json.recon_parameters.phc_flag;       % 1=yes, 0=no
gridding_flag = json.recon_parameters.gridding_flag;  % 1=yes, 0=no
cfc_flag      = json.recon_parameters.cfc_flag;       % 1=yes, 0=no
sfc_flag      = json.recon_parameters.sfc_flag;       % 1=yes, 0=no
gnc_flag      = json.recon_parameters.gnc_flag;       % 1=yes, 0=no

%% Read a .cfl file
%--------------------------------------------------------------------------
% img (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
img_filename = sprintf('img_maxgirf_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_%s_i%d_l%4.2f', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, slice_type, maxiter, lambda);
cfl_file = fullfile(output_path1, img_filename);
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img1 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file2);
fid = fopen(json_file2);
json_txt = fread(fid, [1 inf], 'char=>char');
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
output_path2 = strrep(json.output_path, '/', '\');
% Hack
output_path2 = strrep(output_path2, 'D:', 'E:\projects_lenovo_20250319');

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
Lmax          = json.recon_parameters.Lmax;           % maximum rank of the SVD approximation of a higher-order encoding matrix
L             = json.recon_parameters.L;              % rank of the SVD approximation of a higher-order encoding matrix
lambda        = json.recon_parameters.lambda;         % l2 regularization parameter
tol           = json.recon_parameters.tol;            % PCG tolerance
maxiter       = json.recon_parameters.maxiter;        % PCG maximum iteration 
slice_type    = json.recon_parameters.slice_type;     % type of an excitation slice: "curved" vs "flat"
phc_flag      = json.recon_parameters.phc_flag;       % 1=yes, 0=no
gridding_flag = json.recon_parameters.gridding_flag;  % 1=yes, 0=no
cfc_flag      = json.recon_parameters.cfc_flag;       % 1=yes, 0=no
sfc_flag      = json.recon_parameters.sfc_flag;       % 1=yes, 0=no
gnc_flag      = json.recon_parameters.gnc_flag;       % 1=yes, 0=no

%% Read a .cfl file
%--------------------------------------------------------------------------
% img (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
img_filename = sprintf('img_maxgirf_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_%s_i%d_l%4.2f', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, slice_type, maxiter, lambda);
cfl_file = fullfile(output_path2, img_filename);
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img2 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file3);
fid = fopen(json_file3);
json_txt = fread(fid, [1 inf], 'char=>char');
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
output_path3 = strrep(json.output_path, '/', '\');
% Hack
output_path3 = strrep(output_path3, 'D:', 'E:\projects_lenovo_20250319');

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
Lmax          = json.recon_parameters.Lmax;           % maximum rank of the SVD approximation of a higher-order encoding matrix
L             = json.recon_parameters.L;              % rank of the SVD approximation of a higher-order encoding matrix
lambda        = json.recon_parameters.lambda;         % l2 regularization parameter
tol           = json.recon_parameters.tol;            % PCG tolerance
maxiter       = json.recon_parameters.maxiter;        % PCG maximum iteration 
slice_type    = json.recon_parameters.slice_type;     % type of an excitation slice: "curved" vs "flat"
phc_flag      = json.recon_parameters.phc_flag;       % 1=yes, 0=no
gridding_flag = json.recon_parameters.gridding_flag;  % 1=yes, 0=no
cfc_flag      = json.recon_parameters.cfc_flag;       % 1=yes, 0=no
sfc_flag      = json.recon_parameters.sfc_flag;       % 1=yes, 0=no
gnc_flag      = json.recon_parameters.gnc_flag;       % 1=yes, 0=no

%% Read a .cfl file
%--------------------------------------------------------------------------
% img (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
img_filename = sprintf('img_maxgirf_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_%s_i%d_l%4.2f', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, slice_type, maxiter, lambda);
cfl_file = fullfile(output_path3, img_filename);
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img3 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file4);
fid = fopen(json_file4);
json_txt = fread(fid, [1 inf], 'char=>char');
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
output_path4 = strrep(json.output_path, '/', '\');
% Hack
output_path4 = strrep(output_path4, 'D:', 'E:\projects_lenovo_20250319');

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
Lmax          = json.recon_parameters.Lmax;           % maximum rank of the SVD approximation of a higher-order encoding matrix
L             = json.recon_parameters.L;              % rank of the SVD approximation of a higher-order encoding matrix
lambda        = json.recon_parameters.lambda;         % l2 regularization parameter
tol           = json.recon_parameters.tol;            % PCG tolerance
maxiter       = json.recon_parameters.maxiter;        % PCG maximum iteration 
slice_type    = json.recon_parameters.slice_type;     % type of an excitation slice: "curved" vs "flat"
phc_flag      = json.recon_parameters.phc_flag;       % 1=yes, 0=no
gridding_flag = json.recon_parameters.gridding_flag;  % 1=yes, 0=no
cfc_flag      = json.recon_parameters.cfc_flag;       % 1=yes, 0=no
sfc_flag      = json.recon_parameters.sfc_flag;       % 1=yes, 0=no
gnc_flag      = json.recon_parameters.gnc_flag;       % 1=yes, 0=no

%% Read a .cfl file
%--------------------------------------------------------------------------
% img (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
img_filename = sprintf('img_maxgirf_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_%s_i%d_l%4.2f', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, slice_type, maxiter, lambda);
cfl_file = fullfile(output_path4, img_filename);
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img4 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% x_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, sprintf('x_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
x = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% y_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, sprintf('y_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
y = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% z_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, sprintf('z_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
z = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dx (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, sprintf('dx_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
dx = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dy (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, sprintf('dy_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
dy = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dz (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, sprintf('dz_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
dz = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% read_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, 'read_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
read_sign = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% phase_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path4, 'phase_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
phase_sign = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Set variables
[Nkx,Nky,Nkz] = size(img1);

[Nx,Ny,Nz] = size(img1_dicom);

N = Nx * Ny * Nz;

%% Adjust the image size of a custom reconstruction
idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(Nkx/2) + 1;
idx2_range = (-floor(Ny/2):ceil(Ny/2)-1).' + floor(Nky/2) + 1;
idx3_range = (-floor(Nz/2):ceil(Nz/2)-1).' + floor(Nkz/2) + 1;

img1 = img1(idx1_range, idx2_range, idx3_range);
img2 = img2(idx1_range, idx2_range, idx3_range);
img3 = img3(idx1_range, idx2_range, idx3_range);
img4 = img4(idx1_range, idx2_range, idx3_range);

x = x(idx1_range, idx2_range, idx3_range);
y = y(idx1_range, idx2_range, idx3_range);
z = z(idx1_range, idx2_range, idx3_range);

dx = dx(idx1_range, idx2_range, idx3_range);
dy = dy(idx1_range, idx2_range, idx3_range);
dz = dz(idx1_range, idx2_range, idx3_range);

%% Flip variables
if read_sign == -1
    img1 = flip(img1,1);
    img2 = flip(img2,1);
    img3 = flip(img3,1);
    img4 = flip(img4,1);

    x = flip(x,1);
    y = flip(y,1);
    z = flip(z,1);

    dx = flip(dx,1);
    dy = flip(dy,1);
    dz = flip(dz,1);
end

if phase_sign == -1
    img1 = flip(img1,2);
    img2 = flip(img2,2);
    img3 = flip(img3,2);
    img4 = flip(img4,2);

    x = flip(x,2);
    y = flip(y,2);
    z = flip(z,2);

    dx = flip(dx,2);
    dy = flip(dy,2);
    dz = flip(dz,2);
end

%% Zeropad in k-space
zpad_factor = 4;

slice_number = 32;

ksp1 = sqrt(Nx * Ny) * fftshift(ifft2(ifftshift(img1(:,:,slice_number))));
ksp2 = sqrt(Nx * Ny) * fftshift(ifft2(ifftshift(img2(:,:,slice_number))));
ksp3 = sqrt(Nx * Ny) * fftshift(ifft2(ifftshift(img3(:,:,slice_number))));
ksp4 = sqrt(Nx * Ny) * fftshift(ifft2(ifftshift(img4(:,:,slice_number))));

ksp1_dicom = sqrt(Nx * Ny) * fftshift(ifft2(ifftshift(img1_dicom(:,:,slice_number))));
ksp2_dicom = sqrt(Nx * Ny) * fftshift(ifft2(ifftshift(img2_dicom(:,:,slice_number))));

idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(zpad_factor * Nx / 2) + 1;
idx2_range = (-floor(Ny/2):ceil(Ny/2)-1).' + floor(zpad_factor * Ny / 2) + 1;

ksp1_zpad = complex(zeros([Nx Ny] * zpad_factor, 'single'));
ksp2_zpad = complex(zeros([Nx Ny] * zpad_factor, 'single'));
ksp3_zpad = complex(zeros([Nx Ny] * zpad_factor, 'single'));
ksp4_zpad = complex(zeros([Nx Ny] * zpad_factor, 'single'));

ksp1_dicom_zpad = complex(zeros([Nx Ny] * zpad_factor, 'single'));
ksp2_dicom_zpad = complex(zeros([Nx Ny] * zpad_factor, 'single'));

ksp1_zpad(idx1_range,idx2_range) = ksp1;
ksp2_zpad(idx1_range,idx2_range) = ksp2;
ksp3_zpad(idx1_range,idx2_range) = ksp3;
ksp4_zpad(idx1_range,idx2_range) = ksp4;

ksp1_dicom_zpad(idx1_range,idx2_range) = ksp1_dicom;
ksp2_dicom_zpad(idx1_range,idx2_range) = ksp2_dicom;

%% Calculate sinc-interpolated images
img1_interp = 1 / sqrt(Nx * Ny) * fftshift(fft2(ifftshift(ksp1_zpad)));
img2_interp = 1 / sqrt(Nx * Ny) * fftshift(fft2(ifftshift(ksp2_zpad)));
img3_interp = 1 / sqrt(Nx * Ny) * fftshift(fft2(ifftshift(ksp3_zpad)));
img4_interp = 1 / sqrt(Nx * Ny) * fftshift(fft2(ifftshift(ksp4_zpad)));

img1_dicom_interp = 1 / sqrt(Nx * Ny) * fftshift(fft2(ifftshift(ksp1_dicom_zpad)));
img2_dicom_interp = 1 / sqrt(Nx * Ny) * fftshift(fft2(ifftshift(ksp2_dicom_zpad)));

%% Scale images
c1 = floor(2*Nx/2) + 1;
c2 = floor(2*Ny/2) + 1;

c1 = 250;
c2 = 210;

scale_factor1 = abs(img1_interp(c1,c2));
scale_factor2 = abs(img2_interp(c1,c2));
scale_factor3 = abs(img3_interp(c1,c2));
scale_factor4 = abs(img4_interp(c1,c2));

scale_factor_dicom1 = abs(img1_dicom_interp(c1,c2));
scale_factor_dicom2 = abs(img2_dicom_interp(c1,c2));

img1_scaled = img1_interp / scale_factor1;
img2_scaled = img2_interp / scale_factor2;
img3_scaled = img3_interp / scale_factor3;
img4_scaled = img4_interp / scale_factor4;

img1_dicom_scaled = img1_dicom_interp / scale_factor_dicom1;
img2_dicom_scaled = img2_dicom_interp / scale_factor_dicom2;

%% Display images
baby_blue = [193 220 243] / 255;
blue      = [0   173 236] / 255;
orange    = [239 173 127] / 255;
green     = [205 235 188] / 255;
yellow    = [253 234 155] / 255;

orange_siemens = [236 99 0] / 255;
green_siemens = [3 153 153] / 255;

red_color = [201 37 31] / 255;
blue_color = [86 120 191] / 255;

color_order{1} = '#1f77b4';
color_order{2} = '#ff7f0e';
color_order{3} = '#2ca02c';
color_order{4} = '#d62728';
color_order{5} = '#9467bd';
color_order{6} = '#8c564b';
color_order{7} = '#e377c2';
color_order{8} = '#7f7f7f';
color_order{9} = '#bcbd22';
color_order{10} = '#17becf';

color_order_rgb = hex2rgb(color_order);

cmap = flip(brewermap([],"RdBu"),1);

FontSize = 12;

climits1 = [0 2.0];
climits2 = [0 1.9];

idx1_range_zoom = (300:800).';
idx2_range_zoom = (550:820).';

N1_zoom = length(idx1_range_zoom);
N2_zoom = length(idx2_range_zoom);

figure('Color', 'w', 'Position', [3 5 950 987]);

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/0
%--------------------------------------------------------------------------
ax1 = subplot(4,4,1);
hold on;
imagesc(ax1, abs(img1_scaled.'));

% left, vertical (up-down)
plot(ax1, [idx1_range_zoom(1) idx1_range_zoom(1)], [idx2_range_zoom(1) idx2_range_zoom(end)], 'Color', red_color, 'LineWidth', 1.5);

% right, vertical (up-down)
plot(ax1, [idx1_range_zoom(end) idx1_range_zoom(end)], [idx2_range_zoom(1) idx2_range_zoom(end)], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax1, [idx1_range_zoom(1) idx1_range_zoom(end)], [idx2_range_zoom(1) idx2_range_zoom(1)], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax1, [idx1_range_zoom(1) idx1_range_zoom(end)], [idx2_range_zoom(end) idx2_range_zoom(end)], 'Color', red_color, 'LineWidth', 1.5);

axis image ij off;
colormap(ax1, gray(256));
clim(ax1, climits1);
title(ax1, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax1, sprintf('Gridding/PHC'), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
text(ax1, 2, 0, '(A)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax1, 0, N2 * zpad_factor / 2, sprintf('Slice at z = %3.1f mm', z(1,1,slice_number) * 1e3), 'Rotation', 90, 'FontSize', FontSize, 'FontWeight', 'normal', 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(ax1, 0, N2 * zpad_factor / 2, sprintf('PE direction (A >> P)'), 'FontSize', FontSize, 'Rotation', -90, 'Interpreter', 'tex', 'FontWeight', 'normal', 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Magnitude: Vendor reconstruction Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/0
%--------------------------------------------------------------------------
ax2 = subplot(4,4,2);
imagesc(ax2, abs(img1_dicom_scaled.'));
axis image off;
colormap(ax2, gray(256));
clim(ax2, climits2);
title(ax2, 'Traditional recon', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax2, sprintf('Gridding/PHC'), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
text(ax2, 2, 0, '(B)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% title
%--------------------------------------------------------------------------
text(ax2, N2 * zpad_factor + 50, -490, sprintf('3D GRE-EPI: axial, 0.78 x 0.78 x 0.90 mm^3, 176 slices, R = 1, no PF, ETL = 61, 1 NSA'), 'Color', blue, 'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

text(ax2, N2 * zpad_factor + 60, -260, {sprintf('Gridding/PHC: gridding for ramp sampling/odd-even echo phase correction'), ...
                      sprintf('{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC} vs {\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC(3D)}: concomitant field correction/gradient nonlinearity correction {\\color[rgb]{%f %f %f}during recon} vs {\\color[rgb]{%f %f %f}after recon}', ...
                      orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), ...
                      orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3))}, ...
                      'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Create line (top)
annotation(gcf, 'line', [0.0612 0.9590], [0.8720 0.8720], 'LineWidth', 2);

% Create line (bottom)
annotation(gcf, 'line', [0.0612 0.9590], [0.8250 0.8250], 'LineWidth', 2);

% Create arrow
annotation(gcf, 'arrow', [0.243157894736842 0.233684210526316],...
    [0.729483282674773 0.718338399189464], 'Color', [1 0 0], 'HeadWidth', 6,...
    'HeadStyle', 'plain', 'HeadLength', 6, 'LineWidth', 2);

% Create arrow
annotation(gcf, 'arrow', [0.455789473684209 0.446315789473683],...
    [0.729483282674773 0.718338399189464], 'Color', [1 0 0], 'HeadWidth', 6,...
    'HeadStyle', 'plain', 'HeadLength', 6, 'LineWidth', 2);

% Create arrow (green)
annotation(gcf, 'arrow', [0.587368421052629 0.577894736842103],...
    [0.124620060790277 0.113475177304968], 'Color', color_order{3}, 'HeadWidth', 6,...
    'HeadStyle', 'plain', 'HeadLength', 6, 'LineWidth', 2);

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/0
%--------------------------------------------------------------------------
ax3 = subplot(4,4,3);
imagesc(ax3, abs(img3_scaled.'));
axis image off;
colormap(ax3, gray(256));
clim(ax3, climits1);
title(ax3, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax3, sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}CFC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, 2, 0, '(C)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Displacement field along the x-axis
%--------------------------------------------------------------------------
ax4 = subplot(4,4,4);
hold on;
imagesc(ax4, dx(:,:,slice_number).' * 1e3);
contour(ax4, dx(:,:,slice_number).' * 1e3, (-3:0.5:3).', 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'k');
axis image ij off;
colormap(ax4, cmap);
clim(ax4, [-2.2 2.2]);
title(ax4, 'Displacement field', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax4, {'along the x-axis (RO direction)'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax4, N1 / 2 - 10, 15, {'x-axis $$\longrightarrow$$'}, 'Rotation', 0, 'Color', 'k', 'Interpreter', 'latex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(ax4, 2, 0, '(D)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
hc4 = colorbar;
set(hc4, 'Position', [0.9274 0.5836 0.0126 0.1713], 'FontSize', FontSize);
hTitle4 = title(hc4, '[mm]', 'FontSize', FontSize, 'Position', [11.9887 131.0048 0]);

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/1
%--------------------------------------------------------------------------
ax5 = subplot(4,4,5);
imagesc(ax5, abs(img2_scaled.'));
axis image off;
colormap(ax5, gray(256));
clim(ax5, climits1);
title(ax5, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax5, sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}GNC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax5, 2, 0, '(E)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Magnitude: Vendor reconstruction Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/1
%--------------------------------------------------------------------------
ax6 = subplot(4,4,6);
imagesc(ax6, abs(img2_dicom_scaled.'));
axis image off;
colormap(ax6, gray(256));
clim(ax6, climits2);
title(ax6, 'Traditional recon', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax6, sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}GNC}', green_siemens(1,1), green_siemens(1,2), green_siemens(1,3)), 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax6, 2, 0, '(F)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/1
%--------------------------------------------------------------------------
ax7 = subplot(4,4,7);
imagesc(ax7, abs(img4_scaled.'));
axis image off;
colormap(ax7, gray(256));
clim(ax7, climits1);
title(ax7, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax7, sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax7, 2, 0, '(G)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Displacement field along the y-axis
%--------------------------------------------------------------------------
ax8 = subplot(4,4,8);
hold on;
imagesc(ax8, dy(:,:,slice_number).' * 1e3);
contour(ax8, dy(:,:,slice_number).' * 1e3, (-3:0.5:3).', 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'k');
axis image ij off;
colormap(ax8, cmap);
clim(ax8, [-2.2 2.2]);
title(ax8, 'Displacement field', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax8, {'along the y-axis (PE direction)'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax8, 0, N2/2-12, {'y-axis $$\longrightarrow$$'}, 'Rotation', 90, 'Color', 'k', 'Interpreter', 'latex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(ax8, 2, 0, '(H)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
hc8 = colorbar;
set(hc8, 'Position', [0.9274 0.3374 0.0126 0.1713], 'FontSize', FontSize);
hTitle8 = title(hc8, '[mm]', 'FontSize', FontSize, 'Position', [11.9887 131.0048 0]);

%--------------------------------------------------------------------------
% Zoom: Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/0
%--------------------------------------------------------------------------
ax9 = subplot(4,4,9);
imagesc(ax9, abs(img1_scaled(idx1_range_zoom,idx2_range_zoom).'));
axis image off;
colormap(ax9, gray(256));
clim(ax9, climits1 * 0.8);
text(ax9, 2, 0, '(I)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax9, N1_zoom/2, 0, 'Magnified inset from (A)', 'FontSize', 12, 'FontWeight', 'normal', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Zoom: Vendor reconstruction Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/0
%--------------------------------------------------------------------------
ax10 = subplot(4,4,10);
imagesc(ax10, abs(img1_dicom_scaled(idx1_range_zoom,idx2_range_zoom).'));
axis image off;
colormap(ax10, gray(256));
clim(ax10, climits2 * 0.8);
text(ax10, 2, 0, '(J)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax10, N1_zoom/2, 0, 'Magnified inset from (B)', 'FontSize', 12, 'FontWeight', 'normal', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Zoom: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/0
%--------------------------------------------------------------------------
ax11 = subplot(4,4,11);
imagesc(ax11, abs(img3_scaled(idx1_range_zoom,idx2_range_zoom).'));
axis image off;
colormap(ax11, gray(256));
clim(ax11, climits1 * 0.8);
text(ax11, 2, 0, '(K)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax11, N1_zoom/2, 0, 'Magnified inset from (C)', 'FontSize', 12, 'FontWeight', 'normal', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Zoom: Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/1
%--------------------------------------------------------------------------
ax12 = subplot(4,4,12);
imagesc(ax12, abs(img2_scaled(idx1_range_zoom,idx2_range_zoom).'));
axis image off;
colormap(ax12, gray(256));
clim(ax12, climits1 * 0.8);
text(ax12, 2, 0, '(L)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax12, N1_zoom/2, 0, 'Magnified inset from (E)', 'FontSize', 12, 'FontWeight', 'normal', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Zoom: Vendor reconstruction Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/1
%--------------------------------------------------------------------------
ax13 = subplot(4,4,13);
imagesc(ax13, abs(img2_dicom_scaled(idx1_range_zoom,idx2_range_zoom).'));
axis image off;
colormap(ax13, gray(256));
clim(ax13, climits2 * 0.8);
text(ax13, 2, 0, '(M)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax13, N1_zoom/2, 0, 'Magnified inset from (F)', 'FontSize', 12, 'FontWeight', 'normal', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Zoom: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/1
%--------------------------------------------------------------------------
ax14 = subplot(4,4,14);
imagesc(ax14, abs(img4_scaled(idx1_range_zoom,idx2_range_zoom).'));
axis image off;
colormap(ax14, gray(256));
clim(ax14, climits1 * 0.8);
text(ax14, 2, 0, '(N)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax14, N1_zoom/2, 0, 'Magnified inset from (G)', 'FontSize', 12, 'FontWeight', 'normal', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

set(ax1, 'Position', [0.0768 0.6302-0.1 0.2098 0.2948]);
set(ax2, 'Position', [0.2884 0.6302-0.1 0.2098 0.2948]); %+2116
set(ax3, 'Position', [0.5000 0.6302-0.1 0.2098 0.2948]);
set(ax4, 'Position', [0.7116 0.6302-0.1 0.2098 0.2948]);

set(ax5, 'Position', [0.0768 0.3836-0.1 0.2098 0.2948]);
set(ax6, 'Position', [0.2884 0.3836-0.1 0.2098 0.2948]);
set(ax7, 'Position', [0.5000 0.3836-0.1 0.2098 0.2948]);
set(ax8, 'Position', [0.7116 0.3836-0.1 0.2098 0.2948]);

set(ax9 , 'Position', [0.0763 0.2249-0.1 0.2790 0.2570]);
set(ax10, 'Position', [0.3595 0.2249-0.1 0.2790 0.2570]);
set(ax11, 'Position', [0.6427 0.2249-0.1 0.2790 0.2570]);

set(ax12, 'Position', [0.0763 0.0755-0.1 0.2790 0.2570]);
set(ax13, 'Position', [0.3595 0.0755-0.1 0.2790 0.2570]);
set(ax14, 'Position', [0.6427 0.0755-0.1 0.2790 0.2570]);

export_fig(sprintf('figure7_slc%d', slice_number), '-r300', '-tif', '-c[280, 100, 80, 170]'); % [top,right,bottom,left]
%close gcf;

