% demo_batch_acr_phantom_20250110.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/13/2025, Last modified: 01/13/2025

% Body-6 coil

%% Clean slate
close all; clear all; clc;

%% Set source directories
if ispc
    package_path   = 'D:\cartesian_maxgirf_epi_3d';
    ismrmrd_path   = 'D:\ismrmrd';
    grad_file_path = 'D:\cartesian_maxgirf_epi_3d\GradientCoils';
else
    package_path = '/server/sdata/nlee/cartesian_maxgirf_epi_3d';
    ismrmrd_path = '/server/sdata/nlee/ismrmrd';
    grad_file_path = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/GradientCoils';
end

%% Add source directories to search path
addpath(genpath(package_path));
addpath(genpath(ismrmrd_path));

%% Define a list of .json files
if ispc
    json_files{1} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding0_phc0_cfc0_sfc0_gnc0.json';
    json_files{2} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc0_cfc0_sfc0_gnc0.json';
    json_files{3} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc0_sfc0_gnc0.json';
    json_files{4} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc0_gnc0.json';
    json_files{5} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc1_gnc0.json';
    json_files{6} = 'D:\cartesian_maxgirf_epi_3d\data\acr_phantom_20250110\meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc1_gnc1.json';
else
    json_files{1} = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/data/acr_phantom_20250110/meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding0_phc0_cfc0_sfc0_gnc0_server.json';
    json_files{2} = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/data/acr_phantom_20250110/meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc0_cfc0_sfc0_gnc0_server.json';
    json_files{3} = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/data/acr_phantom_20250110/meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc0_sfc0_gnc0_server.json';
    json_files{4} = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/data/acr_phantom_20250110/meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc0_gnc0_server.json';
    json_files{5} = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/data/acr_phantom_20250110/meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc1_gnc0_server.json';
    json_files{6} = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/data/acr_phantom_20250110/meas_MID00261_FID07209_ep3d_tra_RL_avg1_etl101_gridding1_phc1_cfc1_sfc1_gnc1_server.json';
end

%% Calculate the number of json files
nr_json_files = length(json_files);  

%% Process per json file
for json_number = 5%1:nr_json_files

    %% Define the name of a .json file
    json_file = json_files{json_number};

    %% Calculate voxel coordinates
    demo_step1_cartesian_maxgirf_3d_calculate_voxel_coordinates;

    %% Calculate a fieldmap
    demo_calculate_fieldmap_gre_field_mapping;

    %% Prepare "imaging" k-space data
    demo_step3_cartesian_maxgirf_3d_prepare_ksp_imaging;

    return

    %% Estimate CSMs
    demo_step4_cartesian_maxgirf_3d_estimate_csm;

    %% Cartesian MaxGIRF reconstruction
    demo_step5_cartesian_maxgirf_3d_recon;
end
