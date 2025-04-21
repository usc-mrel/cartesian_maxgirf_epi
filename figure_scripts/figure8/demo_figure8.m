% demo_figure8.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 03/04/2025, Last modified: 03/11/2025

%% Clean slate
close all; clear all; clc;

%% Start a stopwatch timer
start_time = tic;

%% Set source directories
package_path = 'D:\cartesian_maxgirf_epi_3d';

%% Add source directories to search path
addpath(genpath(package_path));
addpath D:\cartesian_maxgirf_epi_3d\figure_scripts\figure_vol1109\tensor_denoising\Tensor-MP-PCA;

%% Define the full path of an output directory
output_path1{1} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00107_FID11278_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan1_gridding1_phc1_cfc0_sfc0_gnc0';
output_path1{2} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00109_FID11280_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan2_gridding1_phc1_cfc0_sfc0_gnc0';
output_path1{3} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00111_FID11282_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan3_gridding1_phc1_cfc0_sfc0_gnc0';

output_path2{1} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00107_FID11278_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan1_gridding1_phc1_cfc1_sfc0_gnc0';
output_path2{2} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00109_FID11280_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan2_gridding1_phc1_cfc1_sfc0_gnc0';
output_path2{3} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00111_FID11282_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan3_gridding1_phc1_cfc1_sfc0_gnc0';

output_path3{1} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00107_FID11278_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan1_gridding1_phc1_cfc1_sfc0_gnc1';
output_path3{2} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00109_FID11280_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan2_gridding1_phc1_cfc1_sfc0_gnc1';
output_path3{3} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00111_FID11282_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan3_gridding1_phc1_cfc1_sfc0_gnc1';

dicom_path1{1} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\dicom\S33_ep3d_tra_AP_highres_TR186_TE80_etl61_0.8mm_scan1'; % ND
dicom_path1{2} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\dicom\S34_ep3d_tra_AP_highres_TR186_TE80_etl61_0.8mm_scan2'; % ND
dicom_path1{3} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\dicom\S35_ep3d_tra_AP_highres_TR186_TE80_etl61_0.8mm_scan3'; % ND

dicom_path2{1} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\dicom\S38_ep3d_tra_AP_highres_TR186_TE80_etl61_0.8mm_scan1_S33_DIS3D'; % DIS3D
dicom_path2{2} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\dicom\S39_ep3d_tra_AP_highres_TR186_TE80_etl61_0.8mm_scan2_S34_DIS3D'; % DIS3D
dicom_path2{3} = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\dicom\S40_ep3d_tra_AP_highres_TR186_TE80_etl61_0.8mm_scan3_S35_DIS3D'; % DIS3D

%E:\projects_lenovo_20250319\cartesian_maxgirf_epi_2d\data

output_path1{1} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00107_FID11278_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan1_gridding1_phc1_cfc0_sfc0_gnc0';
output_path1{2} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00109_FID11280_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan2_gridding1_phc1_cfc0_sfc0_gnc0';
output_path1{3} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00111_FID11282_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan3_gridding1_phc1_cfc0_sfc0_gnc0';

output_path2{1} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00107_FID11278_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan1_gridding1_phc1_cfc1_sfc0_gnc0';
output_path2{2} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00109_FID11280_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan2_gridding1_phc1_cfc1_sfc0_gnc0';
output_path2{3} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00111_FID11282_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan3_gridding1_phc1_cfc1_sfc0_gnc0';

output_path3{1} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00107_FID11278_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan1_gridding1_phc1_cfc1_sfc0_gnc1';
output_path3{2} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00109_FID11280_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan2_gridding1_phc1_cfc1_sfc0_gnc1';
output_path3{3} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00111_FID11282_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan3_gridding1_phc1_cfc1_sfc0_gnc1';

dicom_path1{1} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\dicom\S33_ep3d_tra_AP_highres_TR186_TE80_etl61_0.8mm_scan1'; % ND
dicom_path1{2} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\dicom\S34_ep3d_tra_AP_highres_TR186_TE80_etl61_0.8mm_scan2'; % ND
dicom_path1{3} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\dicom\S35_ep3d_tra_AP_highres_TR186_TE80_etl61_0.8mm_scan3'; % ND

dicom_path2{1} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\dicom\S38_ep3d_tra_AP_highres_TR186_TE80_etl61_0.8mm_scan1_S33_DIS3D'; % DIS3D
dicom_path2{2} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\dicom\S39_ep3d_tra_AP_highres_TR186_TE80_etl61_0.8mm_scan2_S34_DIS3D'; % DIS3D
dicom_path2{3} = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi_3d\data\vol1109_20250302\dicom\S40_ep3d_tra_AP_highres_TR186_TE80_etl61_0.8mm_scan3_S35_DIS3D'; % DIS3D

%% Calculate the number of averages
nr_averages = length(dicom_path1);

%% Get directory information
dir_info = dir(fullfile(dicom_path1{1}, '*IMA'));
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

%% Read DICOM files
img1_dicom = zeros(N1, N2, nr_files, nr_averages, 'single');
x_dicom1 = zeros(N1, N2, nr_files, 'single');
y_dicom1 = zeros(N1, N2, nr_files, 'single');
z_dicom1 = zeros(N1, N2, nr_files, 'single');

for avg = 1:nr_averages

    %% Get directory information
    dir_info = dir(fullfile(dicom_path1{avg}, '*IMA'));
    nr_files = length(dir_info);

    for idx = 1:nr_files

        %% Set the fullpath of a DICOM file
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

        img1_dicom(:,:,idx,avg) = RescaleSlope * double(dicomread(dicom_info)).' + RescaleIntercept; % transpose it!

        %% Parse the DICOM header
        %------------------------------------------------------------------
        % Patient Position Attribute
        % Patient position descriptor relative to the equipment
        %------------------------------------------------------------------
        patient_position = dicom_info.PatientPosition;

        %------------------------------------------------------------------
        % Slice Thickness Attribute
        % Nominal slice thickness, in mm
        %------------------------------------------------------------------
        slice_thickness = dicom_info.SliceThickness; % [mm]

        %------------------------------------------------------------------
        % Image Position (Patient) Attribute
        % The x, y, and z coordinates of the upper left hand corner
        % (center of the first voxel transmitted) of the image, in mm
        %------------------------------------------------------------------
        ipp = dicom_info.ImagePositionPatient; % [mm]

        %------------------------------------------------------------------
        % Image Orientation (Patient) Attribute
        % The direction cosines of the first row and the first column with respect
        % to the patient
        %------------------------------------------------------------------
        iop = dicom_info.ImageOrientationPatient;

        %------------------------------------------------------------------
        % Pixel Spacing Attribute
        % Physical distance in the Patient between the center of each pixel, specified
        % by a numeric pair - adjacent row spacing, adjacent column spacing in mm
        %------------------------------------------------------------------
        pixel_spacing = dicom_info.PixelSpacing; % [mm]

        %------------------------------------------------------------------
        % Number of slices
        %------------------------------------------------------------------
        N3 = 1;

        %------------------------------------------------------------------
        % Slice number
        %------------------------------------------------------------------
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
end

%% Read DICOM files
img2_dicom = zeros(N1, N2, nr_files, nr_averages, 'single');
x_dicom2 = zeros(N1, N2, nr_files, 'single');
y_dicom2 = zeros(N1, N2, nr_files, 'single');
z_dicom2 = zeros(N1, N2, nr_files, 'single');

for avg = 1:nr_averages

    %% Get directory information
    dir_info = dir(fullfile(dicom_path2{avg}, '*IMA'));
    nr_files = length(dir_info);

    for idx = 1:nr_files

        %% Set the fullpath of a DICOM file
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

        img2_dicom(:,:,idx,avg) = RescaleSlope * double(dicomread(dicom_info)).' + RescaleIntercept; % transpose it!

        %% Parse the DICOM header
        %------------------------------------------------------------------
        % Patient Position Attribute
        % Patient position descriptor relative to the equipment
        %------------------------------------------------------------------
        patient_position = dicom_info.PatientPosition;

        %------------------------------------------------------------------
        % Slice Thickness Attribute
        % Nominal slice thickness, in mm
        %------------------------------------------------------------------
        slice_thickness = dicom_info.SliceThickness; % [mm]

        %------------------------------------------------------------------
        % Image Position (Patient) Attribute
        % The x, y, and z coordinates of the upper left hand corner
        % (center of the first voxel transmitted) of the image, in mm
        %------------------------------------------------------------------
        ipp = dicom_info.ImagePositionPatient; % [mm]

        %------------------------------------------------------------------
        % Image Orientation (Patient) Attribute
        % The direction cosines of the first row and the first column with respect
        % to the patient
        %------------------------------------------------------------------
        iop = dicom_info.ImageOrientationPatient;

        %------------------------------------------------------------------
        % Pixel Spacing Attribute
        % Physical distance in the Patient between the center of each pixel, specified
        % by a numeric pair - adjacent row spacing, adjacent column spacing in mm
        %------------------------------------------------------------------
        pixel_spacing = dicom_info.PixelSpacing; % [mm]

        %------------------------------------------------------------------
        % Number of slices
        %------------------------------------------------------------------
        N3 = 1;

        %------------------------------------------------------------------
        % Slice number
        %------------------------------------------------------------------
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
end

%% Set variables
Nkx = 512;
Nky = 306;
Nkz = 176;

[Nx,Ny,Nz,nr_averages] = size(img1_dicom);

%% Read a .cfl file
img1 = complex(zeros(Nkx, Nky, Nkz, nr_averages, 'single'));

for avg = 1:nr_averages
    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path1{avg}, 'img_maxgirf_gridding1_phc1_cfc0_sfc0_gnc0_flat_i6_l0.00');
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    img1(:,:,:,avg) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Read a .cfl file
img2 = complex(zeros(Nkx, Nky, Nkz, nr_averages, 'single'));

for avg = 1:nr_averages
    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path2{avg}, 'img_maxgirf_gridding1_phc1_cfc1_sfc0_gnc0_flat_i6_l0.00');
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    img2(:,:,:,avg) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Read a .cfl file
img3 = complex(zeros(Nkx, Nky, Nkz, nr_averages, 'single'));

for avg = 1:nr_averages
    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path3{avg}, 'img_maxgirf_gridding1_phc1_cfc1_sfc0_gnc1_flat_i6_l0.00');
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    img3(:,:,:,avg) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Read a .cfl file
output_path = 'D:\cartesian_maxgirf_epi_3d\data\vol1109_20250302\meas_MID00111_FID11282_ep3d_tra_AP_highres_TR186_TE80_etl61_0_8mm_scan3_gridding1_phc1_cfc1_sfc0_gnc1';
% Hack
output_path = strrep(output_path, 'D:', 'E:\projects_lenovo_20250319');

slice_type = 'flat';

%--------------------------------------------------------------------------
% x_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('x_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
x = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% y_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('y_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
y = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% z_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('z_shift_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
z = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dx (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('dx_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
dx = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dy (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('dy_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
dy = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% dz (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('dz_%s', slice_type));
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
dz = readcfl(cfl_file);
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

%% Adjust the image size of a custom reconstruction
idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(Nkx/2) + 1;
idx2_range = (-floor(Ny/2):ceil(Ny/2)-1).' + floor(Nky/2) + 1;
idx3_range = (-floor(Nz/2):ceil(Nz/2)-1).' + floor(Nkz/2) + 1;

img1 = img1(idx1_range, idx2_range, idx3_range, :);
img2 = img2(idx1_range, idx2_range, idx3_range, :);
img3 = img3(idx1_range, idx2_range, idx3_range, :);

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

    x = flip(x,2);
    y = flip(y,2);
    z = flip(z,2);

    dx = flip(dx,2);
    dy = flip(dy,2);
    dz = flip(dz,2);
end

%% Read a .cfl file
output_path = 'D:\cartesian_maxgirf_epi_3d\figure_scripts\figure_vol1109';
% Hack
output_path = strrep(output_path, 'D:', 'E:\projects_lenovo_20250319');

%--------------------------------------------------------------------------
% img_denoised (Nx x Ny x Nz x nr_averages)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'img1_denoised');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img1_denoised = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% img_denoised (Nx x Ny x Nz x nr_averages)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'img2_denoised');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img2_denoised = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% img_denoised (Nx x Ny x Nz x nr_averages)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'img3_denoised');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
img3_denoised = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform averaging
img1_avg = mean(img1,4); % Gridding/PHC
img2_avg = mean(img2,4); % Gridding/PHC/CFC
img3_avg = mean(img3,4); % Gridding/PHC/CFC/GNC

img1_denoised_avg = mean(img1_denoised,4); % Gridding/PHC
img2_denoised_avg = mean(img2_denoised,4); % Gridding/PHC/CFC
img3_denoised_avg = mean(img3_denoised,4); % Gridding/PHC/CFC/GNC

img1_dicom_avg = mean(img1_dicom,4); % Gridding/PHC
img2_dicom_avg = mean(img2_dicom,4); % Gridding/PHC/CFC/GNC

%% Scale images
r = 129;
c = 123;
s = 76;

scale_factor1 = abs(img1_avg(r,c,s));
scale_factor2 = abs(img2_avg(r,c,s));
scale_factor3 = abs(img3_avg(r,c,s));

scale_factor1_denoised = abs(img1_denoised_avg(r,c,s));
scale_factor2_denoised = abs(img2_denoised_avg(r,c,s));
scale_factor3_denoised = abs(img3_denoised_avg(r,c,s));

scale_factor1_dicom = abs(img1_dicom_avg(r,c,s));
scale_factor2_dicom = abs(img2_dicom_avg(r,c,s));

img1_avg_scaled = img1_avg / scale_factor1;
img2_avg_scaled = img2_avg / scale_factor2;
img3_avg_scaled = img3_avg / scale_factor3;

img1_denoised_avg_scaled = img1_denoised_avg / scale_factor1_denoised;
img2_denoised_avg_scaled = img2_denoised_avg / scale_factor2_denoised;
img3_denoised_avg_scaled = img3_denoised_avg / scale_factor3_denoised;

img1_dicom_avg_scaled = img1_dicom_avg / scale_factor1_dicom;
img2_dicom_avg_scaled = img2_dicom_avg / scale_factor2_dicom;

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

c1 = floor(Nx/2) + 1;
c2 = floor(Ny/2) + 1;
c3 = floor(Nz/2) + 1;

cmap = flip(brewermap([],"RdBu"),1);

FontSize = 12;

% 107
% 109
axial_slice_number = 108;

idx1_range = (14:235).'; % left-right
idx2_range = (14:256).'; % up-down

N1 = length(idx1_range);
N2 = length(idx2_range);

idx1_range_zoom = (45:118).'- 7; % left-right
idx2_range_zoom = (140:220).'; % up-down

N1_zoom = length(idx1_range_zoom);
N2_zoom = length(idx2_range_zoom);

figure('Color', 'w', 'Position', [3 5 950 987]);

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/0
%--------------------------------------------------------------------------
ax1 = subplot(3,4,1);
hold on;
imagesc(ax1, abs(img1_avg_scaled(idx1_range,idx2_range,axial_slice_number).'));
axis image ij off;
colormap(ax1, gray(256));
clim(ax1, [0 1]);
title(ax1, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax1, sprintf('Gridding/PHC'), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
text(ax1, 2, 0, '(A)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax1, 0, N2 / 2, sprintf('Slice at z = %4.2f mm', z(1,1,axial_slice_number) * 1e3), 'Rotation', 90, 'FontSize', FontSize, 'FontWeight', 'normal', 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(ax1, 0, N2 / 2, sprintf('PE direction (A >> P)'), 'FontSize', FontSize, 'Rotation', -90, 'Interpreter', 'tex', 'FontWeight', 'normal', 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

% left, vertical (up-down)
plot(ax1, [idx1_range_zoom(1) idx1_range_zoom(1)] - idx1_range(1), [idx2_range_zoom(1) idx2_range_zoom(end)], 'Color', red_color, 'LineWidth', 1.5);

% right, vertical (up-down)
plot(ax1, [idx1_range_zoom(end) idx1_range_zoom(end)], [idx2_range_zoom(1) idx2_range_zoom(end)], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax1, [idx1_range_zoom(1) - idx1_range(1) idx1_range_zoom(end)], [idx2_range_zoom(1) idx2_range_zoom(1)], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax1, [idx1_range_zoom(1) - idx1_range(1) idx1_range_zoom(end)], [idx2_range_zoom(end) idx2_range_zoom(end)], 'Color', red_color, 'LineWidth', 1.5);

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/1
%--------------------------------------------------------------------------
ax2 = subplot(3,4,2);
imagesc(ax2, abs(img3_avg_scaled(idx1_range,idx2_range,axial_slice_number).'));
axis image off;
colormap(ax2, gray(256));
clim(ax2, [0 1]);
title(ax2, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax2, sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax2, 2, 0, '(B)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% title
%--------------------------------------------------------------------------
text(ax2, N1, -90, sprintf('3D GRE-EPI: axial, 0.78 x 0.78 x 0.90 mm^3, 176 slices, R = 1, no PF, ETL = 61, 3 NSA'), 'Color', blue, 'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

text(ax2, N1, -45.5, {sprintf('Gridding/PHC: gridding for ramp sampling/odd-even echo phase correction'), ...
                      sprintf('{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC} vs {\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC(3D)}: concomitant field correction/gradient nonlinearity correction {\\color[rgb]{%f %f %f}during recon} vs {\\color[rgb]{%f %f %f}after recon}', ...
                      orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), ...
                      orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3))}, ...
                      'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Create line (top)
annotation(gcf, 'line', [0.0413 0.9383], [0.9393 0.9393], 'LineWidth', 2);

% Create line (bottom)
annotation(gcf, 'line', [0.0413 0.9383], [0.8950 0.8950], 'LineWidth', 2);

%--------------------------------------------------------------------------
% Magnitude: Vendor reconstruction Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/0
%--------------------------------------------------------------------------
ax3 = subplot(3,4,3);
imagesc(ax3, abs(img1_dicom_avg_scaled(idx1_range,idx2_range,axial_slice_number).'));
axis image off;
colormap(ax3, gray(256));
clim(ax3, [0 1]);
title(ax3, 'Traditional reconstruction', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax3, sprintf('Gridding/PHC'), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
text(ax3, 2, 0, '(C)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Magnitude: Vendor reconstruction Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/1
%--------------------------------------------------------------------------
ax4 = subplot(3,4,4);
imagesc(ax4, abs(img2_dicom_avg_scaled(idx1_range,idx2_range,axial_slice_number).'));
axis image off;
colormap(ax4, gray(256));
clim(ax4, [0 1]);
title(ax4, 'Traditional reconstruction', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax4, sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}GNC}', green_siemens(1,1), green_siemens(1,2), green_siemens(1,3)), 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax4, 2, 0, '(D)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/0
%--------------------------------------------------------------------------
ax5 = subplot(3,4,5);
imagesc(ax5, abs(img1_avg_scaled(idx1_range_zoom,idx2_range_zoom,axial_slice_number).'));
axis image off;
colormap(ax5, gray(256));
clim(ax5, [0 1]);
text(ax5, 2, 0, '(E)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax5, 0, N2_zoom / 2, 'Zoomed-in images', 'Rotation', 90, 'FontSize', FontSize, 'FontWeight', 'normal', 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Magnitude: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/1
%--------------------------------------------------------------------------
ax6 = subplot(3,4,6);
imagesc(ax6, abs(img3_avg_scaled(idx1_range_zoom,idx2_range_zoom,axial_slice_number).'));
axis image off;
colormap(ax6, gray(256));
clim(ax6, [0 1]);
text(ax6, 2, 0, '(F)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Create arrow (red)
annotation(gcf, 'arrow', [0.295789473684211 0.309473684210526], ... % location
    [0.477203647416413 0.488348530901728], 'Color', [1 0 0], 'LineWidth', 2,...  arrow direction
    'HeadWidth', 6, 'HeadStyle', 'plain', 'HeadLength', 6);

% Create arrow (red)
annotation(gcf, 'arrow', [0.756842105263158 0.770526315789472], ...
    [0.477203647416413 0.488348530901728], 'Color', [1 0 0], 'LineWidth', 2,...
    'HeadWidth', 6, 'HeadStyle', 'plain', 'HeadLength', 6);

% Create arrow (green)
annotation(gcf, 'arrow', [0.36 0.373684210526315], ...
    [0.392097264437692 0.403242147923007], 'Color', color_order{3}, 'LineWidth', 2,...
    'HeadWidth', 6, 'HeadStyle', 'plain', 'HeadLength', 6);

% Create arrow (green)
annotation(gcf, 'arrow', [0.823157894736837 0.836842105263152], ...
    [0.392097264437692 0.403242147923007], 'Color', color_order{3}, 'LineWidth', 2,...
    'HeadWidth', 6, 'HeadStyle', 'plain', 'HeadLength', 6);

%--------------------------------------------------------------------------
% Magnitude: Vendor reconstruction Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/0
%--------------------------------------------------------------------------
ax7 = subplot(3,4,7);
imagesc(ax7, abs(img1_dicom_avg_scaled(idx1_range_zoom,idx2_range_zoom,axial_slice_number).'));
axis image off;
colormap(ax7, gray(256));
clim(ax7, [0 1]);
text(ax7, 2, 0, '(G)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Magnitude: Vendor reconstruction Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/1
%--------------------------------------------------------------------------
ax8 = subplot(3,4,8);
imagesc(ax8, abs(img2_dicom_avg_scaled(idx1_range_zoom,idx2_range_zoom,axial_slice_number).'));
axis image off;
colormap(ax8, gray(256));
clim(ax8, [0 1]);
text(ax8, 2, 0, '(H)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Phase: Gridding/PHC/CFC/SFC/GNC = 1/1/0/0/0
%--------------------------------------------------------------------------
ax9 = subplot(3,4,9);
imagesc(ax9, angle(img1_avg_scaled(idx1_range,idx2_range,axial_slice_number).') * 180 / pi);
axis image off;
colormap(ax9, hsv(256));
clim(ax9, [-180 180]);
title(ax9, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax9, sprintf('Gridding/PHC'), 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal');
text(ax9, 2, 0, '(I)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Phase: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/1
%--------------------------------------------------------------------------
ax10 = subplot(3,4,10);
imagesc(ax10, angle(img3_avg_scaled(idx1_range,idx2_range,axial_slice_number).') * 180 / pi);
axis image off;
colormap(ax10, hsv(256));
clim(ax10, [-180 180]);
title(ax10, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax10, sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax10, 2, 0, '(J)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
hc10 = colorbar;
set(hc10, 'Location', 'SouthOutside', 'FontSize', FontSize, 'Position', [0.0779 0.0689 0.3411 0.0100]);
hTitle10 = title(hc10, '[deg]', 'FontSize', FontSize, 'Position', [260.2669 -3.3975 0]);

%--------------------------------------------------------------------------
% Displacement field along the x-axis
%--------------------------------------------------------------------------
ax11 = subplot(3,4,11);
hold on;
imagesc(ax11, dx(idx1_range,idx2_range,axial_slice_number).' * 1e3);
contour(ax11, dx(idx1_range,idx2_range,axial_slice_number).' * 1e3, (-3:0.5:3).', 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'k');

% left, vertical (up-down)
plot(ax11, [idx1_range_zoom(1) idx1_range_zoom(1)] - idx1_range(1), [idx2_range_zoom(1) idx2_range_zoom(end)], 'Color', red_color, 'LineWidth', 1.5);

% right, vertical (up-down)
plot(ax11, [idx1_range_zoom(end) idx1_range_zoom(end)], [idx2_range_zoom(1) idx2_range_zoom(end)], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax11, [idx1_range_zoom(1) - idx1_range(1) idx1_range_zoom(end)], [idx2_range_zoom(1) idx2_range_zoom(1)], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax11, [idx1_range_zoom(1) - idx1_range(1) idx1_range_zoom(end)], [idx2_range_zoom(end) idx2_range_zoom(end)], 'Color', red_color, 'LineWidth', 1.5);

axis image ij off;
colormap(ax11, cmap);
clim(ax11, [-2.2 2.2]);
title(ax11, 'Displacement field', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax11, {'along the x-axis (RO direction)'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax11, 2, 0, '(K)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax11, N1 / 2 - 3, 0, {'x-axis $$\longrightarrow$$'}, 'Rotation', 0, 'Color', 'k', 'Interpreter', 'latex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

%--------------------------------------------------------------------------
% Displacement field along the y-axis
%--------------------------------------------------------------------------
ax12 = subplot(3,4,12);
hold on;
imagesc(ax12, dy(idx1_range,idx2_range,axial_slice_number).' * 1e3);
contour(ax12, dy(idx1_range,idx2_range,axial_slice_number).' * 1e3, (-3:0.5:3).', 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'k');

% left, vertical (up-down)
plot(ax12, [idx1_range_zoom(1) idx1_range_zoom(1)] - idx1_range(1), [idx2_range_zoom(1) idx2_range_zoom(end)], 'Color', red_color, 'LineWidth', 1.5);

% right, vertical (up-down)
plot(ax12, [idx1_range_zoom(end) idx1_range_zoom(end)], [idx2_range_zoom(1) idx2_range_zoom(end)], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax12, [idx1_range_zoom(1) - idx1_range(1) idx1_range_zoom(end)], [idx2_range_zoom(1) idx2_range_zoom(1)], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax12, [idx1_range_zoom(1) - idx1_range(1) idx1_range_zoom(end)], [idx2_range_zoom(end) idx2_range_zoom(end)], 'Color', red_color, 'LineWidth', 1.5);

axis image ij off;
colormap(ax12, cmap);
clim(ax12, [-2.2 2.2]);
title(ax12, 'Displacement field', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax12, {'along the y-axis (PE direction)'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax12, 2, 0, '(L)', 'FontSize', 12, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax12, 0, N2 / 2 - 17, {'y-axis $$\longrightarrow$$'}, 'Rotation', 90, 'Color', 'k', 'Interpreter', 'latex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

hc12 = colorbar;
set(hc12, 'Location', 'SouthOutside', 'FontSize', FontSize, 'Position', [0.5411 0.0689 0.3411 0.0100]);
hTitle12 = title(hc12, '[mm]', 'FontSize', FontSize, 'Position', [260.2669 -3.3975 0]);

set(ax1, 'Position', [0.0184+0.01 0.5218 0.2289 0.4184]);
set(ax2, 'Position', [0.2500+0.01 0.5218 0.2289 0.4184]);
set(ax3, 'Position', [0.4816+0.01 0.5218 0.2289 0.4184]);
set(ax4, 'Position', [0.7132+0.01 0.5218 0.2289 0.4184]);

set(ax5, 'Position', [0.0184+0.01 0.2780 0.2289 0.4184]);
set(ax6, 'Position', [0.2500+0.01 0.2780 0.2289 0.4184]);
set(ax7, 'Position', [0.4816+0.01 0.2780 0.2289 0.4184]);
set(ax8, 'Position', [0.7132+0.01 0.2780 0.2289 0.4184]);

set(ax9 , 'Position', [0.0184+0.01 0.0342-0.039 0.2289 0.4184]);
set(ax10, 'Position', [0.2500+0.01 0.0342-0.039 0.2289 0.4184]);
set(ax11, 'Position', [0.4816+0.01 0.0342-0.039 0.2289 0.4184]);
set(ax12, 'Position', [0.7132+0.01 0.0342-0.039 0.2289 0.4184]);

export_fig(sprintf('figure8_slc%d', axial_slice_number), '-r300', '-tif', '-c[80, 120, 140, 20]'); % [top,right,bottom,left]
close gcf;
