% demo_figure6.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/28/2025, Last modified: 03/10/2025

%% Clean slate
close all; clear all; clc;

%% Start a stopwatch timer
start_time = tic;

%% Set source directories
package_path = 'D:\cartesian_maxgirf_epi_3d';

%% Add source directories to search path
addpath(genpath(package_path));

%% Define the full path of an output directory
json_file1 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding0_phc0_cfc0_sfc0_gnc0.json';
json_file2 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc0_cfc0_sfc0_gnc0.json';
json_file3 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc0_sfc0_gnc0.json';
json_file4 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc0_gnc0.json';
json_file5 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc1_gnc0.json';
json_file6 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc1_gnc1.json';

dicom_path1 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\dicom\S7_ep3d_tra_RL_avg1_etl101_ND'; % ND
dicom_path2 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\dicom\S12_ep3d_tra_RL_avg1_etl101_ND_S7_DIS3D'; % DIS3D
%dicom_path2 = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\dicom\S8_ep3d_tra_RL_avg1_etl101'; % DIS2D

%E:\projects_lenovo_20250319\cartesian_maxgirf_epi_2d\data

json_file1 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding0_phc0_cfc0_sfc0_gnc0.json';
json_file2 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc0_cfc0_sfc0_gnc0.json';
json_file3 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc0_sfc0_gnc0.json';
json_file4 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc0_gnc0.json';
json_file5 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc1_gnc0.json';
json_file6 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc1_gnc1.json';

dicom_path1 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\dicom\S7_ep3d_tra_RL_avg1_etl101_ND'; % ND
dicom_path2 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\dicom\S12_ep3d_tra_RL_avg1_etl101_ND_S7_DIS3D'; % DIS3D

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

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file5);
fid = fopen(json_file5);
json_txt = fread(fid, [1 inf], 'char=>char');
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
output_path5 = strrep(json.output_path, '/', '\');
% Hack
output_path5 = strrep(output_path5, 'D:', 'E:\projects_lenovo_20250319');

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
cfl_file = fullfile(output_path5, img_filename);
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img5 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file6);
fid = fopen(json_file6);
json_txt = fread(fid, [1 inf], 'char=>char');
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
output_path6 = strrep(json.output_path, '/', '\');
% Hack
output_path6 = strrep(output_path6, 'D:', 'E:\projects_lenovo_20250319');

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
cfl_file = fullfile(output_path6, img_filename);
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img6 = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% x_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, sprintf('x_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
x = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% y_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, sprintf('y_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
y = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% z_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, sprintf('z_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
z = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dx (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, sprintf('dx_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
dx = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dy (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, sprintf('dy_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
dy = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dz (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, sprintf('dz_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
dz = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% fieldmap (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, sprintf('fieldmap_smooth_gnl%d_%s', gnc_flag, slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
fieldmap = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% read_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, 'read_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
read_sign = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% phase_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path6, 'phase_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
phase_sign = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Set variables
[Nkx,Nky,Nkz] = size(img6);

[Nx,Ny,Nz] = size(fieldmap);

N = Nx * Ny * Nz;

%% Adjust the image size of a custom reconstruction
idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(Nkx/2) + 1;
idx2_range = (-floor(Ny/2):ceil(Ny/2)-1).' + floor(Nky/2) + 1;
idx3_range = (-floor(Nz/2):ceil(Nz/2)-1).' + floor(Nkz/2) + 1;

x = x(idx1_range, idx2_range, idx3_range);
y = y(idx1_range, idx2_range, idx3_range);
z = z(idx1_range, idx2_range, idx3_range);

dx = dx(idx1_range, idx2_range, idx3_range);
dy = dy(idx1_range, idx2_range, idx3_range);
dz = dz(idx1_range, idx2_range, idx3_range);

img1 = img1(idx1_range, idx2_range, idx3_range);
img2 = img2(idx1_range, idx2_range, idx3_range);
img3 = img3(idx1_range, idx2_range, idx3_range);
img4 = img4(idx1_range, idx2_range, idx3_range);
img5 = img5(idx1_range, idx2_range, idx3_range);
img6 = img6(idx1_range, idx2_range, idx3_range);

%% Flip variables
if read_sign == -1
    x = flip(x,1);
    y = flip(y,1);
    z = flip(z,1);

    dx = flip(dx,1);
    dy = flip(dy,1);
    dz = flip(dz,1);

    img1 = flip(img1,1);
    img2 = flip(img2,1);
    img3 = flip(img3,1);
    img4 = flip(img4,1);
    img5 = flip(img5,1);
    img6 = flip(img6,1);

    fieldmap = flip(fieldmap,1);
end

if phase_sign == -1
    x = flip(x,2);
    y = flip(y,2);
    z = flip(z,2);

    dx = flip(dx,2);
    dy = flip(dy,2);
    dz = flip(dz,2);

    img1 = flip(img1,2);
    img2 = flip(img2,2);
    img3 = flip(img3,2);
    img4 = flip(img4,2);
    img5 = flip(img5,2);
    img6 = flip(img6,2);

    fieldmap = flip(fieldmap,2);
end

%% Display coronal images
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

cmap1 = flip(brewermap([],"RdBu"),1);
cmap2 = colorcet('D1');

FontSize = 12;

slice_number = 129;

ylimits = [10 240];

figure('Color', 'w', 'Position', [3 5 950 987]);

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 0/0/0/0/0
%--------------------------------------------------------------------------
ax1 = subplot(3,5,1);
imagesc(ax1, abs(squeeze(img1(slice_number,:,:))));
axis image off;
colormap(ax1, gray(256));
title(ax1, 'Cartesian MaxGIRF', 'FontSize', FontSize);
subtitle(ax1, sprintf('no correction'), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax1, ylimits);
text(ax1, 0 , N2 / 2, sprintf('PE direction (R >> L)'), 'FontSize', FontSize, 'Rotation', -90, 'Interpreter', 'tex', 'FontWeight', 'normal', 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(ax1, 2, ylimits(1), '(A)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax1, 0, N2 / 2, sprintf('Slice at y = %4.2fmm', y(slice_number,1,1) * 1e3), 'FontSize', FontSize, 'Rotation', 90, 'Interpreter', 'tex', 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/0/0/0/0
%--------------------------------------------------------------------------
ax2 = subplot(3,5,2);
imagesc(ax2, abs(squeeze(img2(slice_number,:,:))));
axis image off;
colormap(ax2, gray(256));
title(ax2, 'Cartesian MaxGIRF', 'FontSize', FontSize);
subtitle(ax2, sprintf('Gridding'), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax2, ylimits);
text(ax2, 2, ylimits(1), '(B)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/0
%--------------------------------------------------------------------------
ax3 = subplot(3,5,3);
imagesc(ax3, abs(squeeze(img3(slice_number,:,:))));
axis image off;
colormap(ax3, gray(256));
title(ax3, 'Cartesian MaxGIRF', 'FontSize', FontSize);
subtitle(ax3, sprintf('Gridding/PHC'), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax3, ylimits);
text(ax3, 2, ylimits(1), '(C)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% title
%--------------------------------------------------------------------------
hText1 = text(ax3, Nz / 2, -92, sprintf('3D GRE-EPI: axial, 1.0 mm^3, 160 slices, R = 1, no partial Fourier, ETL = 101, 1 NSA'), 'Color', blue, 'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

hText2 = text(ax3, Nz / 2, -33, {sprintf('Gridding/PHC: gridding for ramp sampling/odd-even echo phase correction'), ...
                        sprintf('{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC} vs {\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC(3D)}: concomitant field correction/gradient nonlinearity correction {\\color[rgb]{%f %f %f}during recon} vs {\\color[rgb]{%f %f %f}after recon}', ...
                        orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), ...
                        orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3)), ...
                        sprintf('{\\color[rgb]{%f %f %f}SFC}: static off-resonance correction {\\color[rgb]{%f %f %f}during recon}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3))}, ...
                        'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Create line (top)
hLine2 = annotation(gcf, 'line', [0.0465 0.9429], [0.9685 0.9685], 'LineWidth', 2);

% Create line (bottom)
hLine1 = annotation(gcf, 'line', [0.0465 0.9429], [0.9045 0.9045], 'LineWidth', 2);

%--------------------------------------------------------------------------
% Magnitude: Vendor reconstruction Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/0
%--------------------------------------------------------------------------
ax4 = subplot(3,5,4);
imagesc(ax4, abs(squeeze(img1_dicom(:,slice_number,:))));
axis image off;
colormap(ax4, gray(256));
title(ax4, 'Traditional recon', 'FontSize', FontSize, 'HorizontalAlignment', 'center');
subtitle(ax4, sprintf('Gridding/PHC'), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax4, ylimits);
text(ax4, 2, ylimits(1), '(D)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
clim([175 1000]);

%--------------------------------------------------------------------------
% Magnitude: Vendor reconstruction Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/1
%--------------------------------------------------------------------------
ax5 = subplot(3,5,5);
imagesc(ax5, abs(squeeze(img2_dicom(:,slice_number,:))));
axis image off;
colormap(ax5, gray(256));
title(ax5, 'Traditional recon', 'FontSize', FontSize, 'HorizontalAlignment', 'center');
subtitle(ax5, sprintf('G./P./{\\color[rgb]{%f %f %f}GNC}', green_siemens(1,1), green_siemens(1,2), green_siemens(1,3)), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax5, ylimits);
text(ax5, 2, ylimits(1), '(E)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
clim([175 1000]);

text(ax5, Nz/2, 7  , 'R', 'FontSize', FontSize, 'Interpreter', 'tex', 'FontWeight', 'bold', 'Color', 'r', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
text(ax5, Nz/2, 240, 'L', 'FontSize', FontSize, 'Interpreter', 'tex', 'FontWeight', 'bold', 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(ax5, 2 , N2 / 2, 'I', 'FontSize', FontSize, 'Interpreter', 'tex', 'FontWeight', 'bold', 'Color', 'r', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
text(ax5, Nz, N2 / 2, 'H', 'FontSize', FontSize, 'Interpreter', 'tex', 'FontWeight', 'bold', 'Color', 'r', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/0
%--------------------------------------------------------------------------
ax6 = subplot(3,5,6);
imagesc(ax6, abs(squeeze(img4(slice_number,:,:))));
axis image off;
colormap(ax6, gray(256));
title(ax6, 'Cartesian MaxGIRF', 'FontSize', FontSize);
subtitle(ax6, sprintf('G./P./{\\color[rgb]{%f %f %f}CFC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax6, ylimits);
text(ax6, 2, ylimits(1), '(F)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Create arrow (red, top)
annotation(gcf, 'arrow', [0.226315789473684 0.213684210526315], ... % location
    [0.521783181357649 0.509646705343919], 'Color', [1 0 0], 'LineWidth', 2,...  arrow direction
    'HeadWidth', 6, 'HeadStyle', 'plain', 'HeadLength', 6);

% Create arrow (red, bottom)
annotation(gcf, 'arrow', [0.226315789473684 0.213684210526315], ... % location
    [0.521783181357649 0.509646705343919] - 0.0993, 'Color', [1 0 0], 'LineWidth', 2,...  arrow direction
    'HeadWidth', 6, 'HeadStyle', 'plain', 'HeadLength', 6);

% Create arrow (green, bottom)
annotation(gcf, 'arrow', [0.0778947368421052 0.068421052631578], ... % location
    [0.328267477203648 0.341460281838348], 'Color', color_order{3}, 'LineWidth', 2,...  arrow direction
    'HeadWidth', 6, 'HeadStyle', 'plain', 'HeadLength', 6);

% Create arrow (green, top)
annotation(gcf, 'arrow', [0.0778947368421052 0.068421052631578], ... % location
    [0.526849037487343 0.540041842122043], 'Color', color_order{3}, 'LineWidth', 2,...  arrow direction
    'HeadWidth', 6, 'HeadStyle', 'plain', 'HeadLength', 6);

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/1/1/1/0
%--------------------------------------------------------------------------
ax7 = subplot(3,5,7);
imagesc(ax7, abs(squeeze(img5(slice_number,:,:))));
axis image off;
colormap(ax7, gray(256));
title(ax7, 'Cartesian MaxGIRF', 'FontSize', FontSize);
subtitle(ax7, sprintf('G./P./{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}SFC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax7, ylimits);
text(ax7, 2, ylimits(1), '(G)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/1/1/1/1
%--------------------------------------------------------------------------
ax8 = subplot(3,5,8);
imagesc(ax8, abs(squeeze(img6(slice_number,:,:))));
axis image off;
colormap(ax8, gray(256));
title(ax8, 'Cartesian MaxGIRF', 'FontSize', FontSize);
subtitle(ax8, sprintf('G./P./{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}SFC}/{\\color[rgb]{%f %f %f}GNC}', ...
              orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax8, ylimits);
text(ax8, 2, ylimits(1), '(H)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Displacement field along the x-axis
%--------------------------------------------------------------------------
ax9 = subplot(3,5,9);
hold on;
imagesc(ax9, squeeze(dx(slice_number,:,:) * 1e3));
contour(ax9, squeeze(dx(slice_number,:,:)) * 1e3, (-3:0.5:3).', 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'k');
axis image ij off;
colormap(ax9, cmap1);
clim(ax9, [-2.2 2.2]);
title(ax9, 'Displacement field', 'FontSize', FontSize);
subtitle(ax9, 'along the x-axis (PE)', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax9, ylimits);
text(ax9, 90, 116, {'x-axis $$\longrightarrow$$'}, 'Rotation', -90, 'Color', 'k', 'Interpreter', 'latex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(ax9, 2, ylimits(1), '(I)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Displacement field along the z-axis
%--------------------------------------------------------------------------
ax10 = subplot(3,5,10);
hold on;
imagesc(ax10, squeeze(dz(slice_number,:,:) * 1e3));
contour(ax10, squeeze(dz(slice_number,:,:)) * 1e3, (-3:0.5:3).', 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'k');
axis image ij off;
colormap(ax10, cmap1);
clim(ax10, [-2.2 2.2]);
title(ax10, 'Displacement field', 'FontSize', FontSize);
subtitle(ax10, 'along the z-axis (SL)', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax10, ylimits);
text(ax10, 90, 116, {'$$\longleftarrow$$ z-axis'}, 'Rotation', 0, 'Color', 'k', 'Interpreter', 'latex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(ax10, 2, ylimits(1), '(J)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

hc10 = colorbar;
set(hc10, 'Position', [0.9495 0.3739-0.04 0.0126 0.2117], 'FontSize', FontSize);
title(hc10, '[mm]', 'FontSize', FontSize, 'Position', [11.9887 161.65 0]);

%--------------------------------------------------------------------------
% Phase: Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/0
%--------------------------------------------------------------------------
ax11 = subplot(3,5,11);
YData = flip(squeeze(x(slice_number,:,1)),1).' * 1e3;
XData = squeeze(z(slice_number,1,:)).' * 1e3;
imagesc(ax11, 'YData', YData, 'XData', XData, 'CData', angle(squeeze(img3(slice_number,:,:))));
axis image;
colormap(ax11, hsv(256));
set(ax11, 'XDir', 'reverse', 'YDir', 'reverse', 'Box', 'On', 'FontSize', FontSize, 'TickLength', [0.0100 0.0250] * 1, 'TickDir', 'out');
title(ax11, 'Cartesian MaxGIRF', 'FontSize', FontSize);
subtitle(ax11, sprintf('Gridding/PHC'), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax11, [YData(ylimits(1)) YData(ylimits(2))]);
hXLabel11 = xlabel(ax11, 'z [mm]', 'Color', 'k', 'FontSize', FontSize, 'Position', [-106.3531  117.4493 -1]);
hYLabel11 = ylabel(ax11, 'x [mm]', 'Color', 'k', 'FontSize', FontSize, 'Position',  [103.4039 2.5999 -1]);
text(ax11, XData(1), YData(ylimits(1)), '(K)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Phase: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/0
%--------------------------------------------------------------------------
ax12 = subplot(3,5,12);
imagesc(ax12, angle(squeeze(img4(slice_number,:,:))));
axis image off;
colormap(ax12, hsv(256));
title(ax12, 'Cartesian MaxGIRF', 'FontSize', FontSize);
subtitle(ax12, sprintf('G./P./{\\color[rgb]{%f %f %f}CFC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax12, ylimits);
text(ax12, 2, ylimits(1), '(L)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Phase: Gridding/PHC/CFC/SFC/GNC = 1/1/1/1/0
%--------------------------------------------------------------------------
ax13 = subplot(3,5,13);
imagesc(ax13, angle(squeeze(img5(slice_number,:,:))));
axis image off;
colormap(ax13, hsv(256));
title(ax13, 'Cartesian MaxGIRF', 'FontSize', FontSize);
subtitle(ax13, sprintf('G./P./{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}SFC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax13, ylimits);
text(ax13, 2, ylimits(1), '(M)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Phase: Gridding/PHC/CFC/SFC/GNC = 1/1/1/1/1
%--------------------------------------------------------------------------
ax14 = subplot(3,5,14);
imagesc(ax14, angle(squeeze(img6(slice_number,:,:))));
axis image off;
colormap(ax14, hsv(256));
title(ax14, 'Cartesian MaxGIRF', 'FontSize', FontSize);
subtitle(ax14, sprintf('G./P./{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}SFC}/{\\color[rgb]{%f %f %f}GNC}', ...
              orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
ylim(ax14, ylimits);
text(ax14, 2, ylimits(1), '(N)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Static off-resonance map
%--------------------------------------------------------------------------
ax15 = subplot(3,5,15);
imagesc(ax15, (squeeze(fieldmap(slice_number,:,:))));
axis image off;
colormap(ax15, cmap2);
clim([-20 20]);
subtitle(ax15, {'Static', 'off-resonance map'}, 'FontSize', FontSize);
ylim(ax15, ylimits);
text(ax15, 2, ylimits(1), '(O)', 'FontSize', FontSize, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

hc15 = colorbar;
set(hc15, 'Position', [0.9495 0.0841-0.04 0.0126 0.2117], 'FontSize', FontSize);
title(hc15, '[Hz]', 'FontSize', FontSize, 'Position', [8.9887 160.9 0]);

set(ax1, 'Position', [0.0579-0.012 0.6292-0.04 0.1789 0.2958]);
set(ax2, 'Position', [0.2379-0.012 0.6292-0.04 0.1789 0.2958]); %+1800
set(ax3, 'Position', [0.4179-0.012 0.6292-0.04 0.1789 0.2958]);
set(ax4, 'Position', [0.5979-0.012 0.6292-0.04 0.1789 0.2958]);
set(ax5, 'Position', [0.7779-0.012 0.6292-0.04 0.1789 0.2958]);

set(ax6 , 'Position', [0.0579-0.012 0.3394-0.04 0.1789 0.2958]);
set(ax7 , 'Position', [0.2379-0.012 0.3394-0.04 0.1789 0.2958]);
set(ax8 , 'Position', [0.4179-0.012 0.3394-0.04 0.1789 0.2958]);
set(ax9 , 'Position', [0.5979-0.012 0.3394-0.04 0.1789 0.2958]);
set(ax10, 'Position', [0.7779-0.012 0.3394-0.04 0.1789 0.2958]);

set(ax11, 'Position', [0.0579-0.012 0.0496-0.04 0.1789 0.2958]);
set(ax12, 'Position', [0.2379-0.012 0.0496-0.04 0.1789 0.2958]);
set(ax13, 'Position', [0.4179-0.012 0.0496-0.04 0.1789 0.2958]);
set(ax14, 'Position', [0.5979-0.012 0.0496-0.04 0.1789 0.2958]);
set(ax15, 'Position', [0.7779-0.012 0.0496-0.04 0.1789 0.2958]);

export_fig(sprintf('figure6_slc%d', slice_number), '-r300', '-tif', '-c[0, 10, 15, 5]'); % [top,right,bottom,left]
close gcf;
