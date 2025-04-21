% batch_acr_phantom_20250209.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 02/09/2024, Last modified: 02/09/2025

% partial Fourier, R=1 vs partial Fourier, R=2

%% Clean slate
close all; clear all; clc;

%% Set source directories
if ispc
    package_path   = 'D:\cartesian_maxgirf_epi_2d';
    ismrmrd_path   = 'D:\ismrmrd';
    grad_file_path = 'D:\cartesian_maxgirf_epi_2d\GradientCoils';
else
    package_path = '/server/sdata/nlee/cartesian_maxgirf_epi_2d';
    ismrmrd_path = '/server/sdata/nlee/ismrmrd';
    grad_file_path = '/server/sdata/nlee/cartesian_maxgirf_epi_2d/GradientCoils';
end

%% Add source directories to search path
addpath(genpath(package_path));
addpath(genpath(ismrmrd_path));

%% Define a list of .json files
if ispc
    % partial Fourier with R = 1
    %json_files{1} = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\meas_MID00123_FID08926_ep2d_se_bw976_tra_AP_pf_R1_gridding1_phc1_cfc0_sfc0_gnc0.json';
    %json_files{2} = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\meas_MID00123_FID08926_ep2d_se_bw976_tra_AP_pf_R1_gridding1_phc1_cfc1_sfc0_gnc0.json';
    %json_files{3} = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\meas_MID00123_FID08926_ep2d_se_bw976_tra_AP_pf_R1_gridding1_phc1_cfc1_sfc0_gnc1.json';
    %json_files{1} = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\meas_MID00123_FID08926_ep2d_se_bw976_tra_AP_pf_R1_gridding1_phc1_cfc1_sfc1_gnc1.json';

    % partial Fourier with R = 2
    json_files{1} = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\meas_MID00125_FID08928_ep2d_se_bw976_tra_AP_pf_R2_gridding1_phc1_cfc1_sfc0_gnc0.json';
    json_files{2} = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\meas_MID00125_FID08928_ep2d_se_bw976_tra_AP_pf_R2_gridding1_phc1_cfc1_sfc0_gnc1.json';
else
    json_files{1} = '/server/sdata/nlee/cartesian_maxgirf_epi_2d/data/acr_phantom_20250209/meas_MID00125_FID08928_ep2d_se_bw976_tra_AP_pf_R2_gridding1_phc1_cfc1_sfc0_gnc0_server.json';
    json_files{2} = '/server/sdata/nlee/cartesian_maxgirf_epi_2d/data/acr_phantom_20250209/meas_MID00125_FID08928_ep2d_se_bw976_tra_AP_pf_R2_gridding1_phc1_cfc1_sfc0_gnc1_server.json';
end

%% Define the name of a .json file for a fieldmap
if ispc
    fieldmap_json_file = 'D:\gre_field_mapping\data\acr_phantom_20250209\meas_MID00152_FID08955_gre_field_mapping_tra_AP_slices56_3mm_TR2360_gnc1.json';
else
    fieldmap_json_file = '/server/sdata/nlee/gre_field_mapping/data/acr_phantom_20250209/meas_MID00262_FID07210_gre_field_mapping_gnl1_server.json';
end

%% Calculate the number of json files
nr_json_files = length(json_files);

%% Process per json file
for json_number = 1:nr_json_files

    %% Define the name of a .json file
    json_file = json_files{json_number};

    %% Calculate voxel coordinates
    %demo_cartesian_maxgirf_2d_calculate_voxel_coordinates;

    %% Calculate a fieldmap
    %demo_calculate_fieldmap_gre_field_mapping_matching_slices;

    %% Prepare "calibration" k-space data
    %demo_cartesian_maxgirf_2d_prepare_ksp_calibration;

    %% Prepare "imaging" k-space data
    %demo_cartesian_maxgirf_2d_prepare_ksp_imaging;

    %% Estimate CSMs
    %demo_cartesian_maxgirf_2d_estimate_csm;

    %% Cartesian MaxGIRF reconstruction
    demo_cartesian_maxgirf_2d_recon;

    %% Cartesian MaxGIRF reconstruction (low resolution)
    %demo_cartesian_maxgirf_2d_recon_lowres;

    %% Phase-constrained, low-rank and subspace model-based MaxGIRF reconstruction
    %demo_cartesian_maxgirf_2d_recon_phase_constrained_lowrank;
end
