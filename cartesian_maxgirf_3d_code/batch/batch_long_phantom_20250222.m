% batch_long_phantom_20250222.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 02/22/2025, Last modified: 02/22/2025

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
    %json_files{1} = 'D:\cartesian_maxgirf_epi_3d\data\long_phantom_20250222\meas_MID00025_FID10270_ep3d_tra_RL_TR217_TE90_etl101_gridding1_phc1_cfc0_sfc0_gnc0.json';
    % seems no ramp sampling is used
    %json_files{2} = 'D:\cartesian_maxgirf_epi_3d\data\long_phantom_20250222\meas_MID00027_FID10272_ep3d_tra_AP_TR79_TE25_etl15_gridding1_phc1_cfc0_sfc0_gnc0.json';
    %json_files{3} = 'D:\cartesian_maxgirf_epi_3d\data\long_phantom_20250222\meas_MID00029_FID10274_ep3d_tra_AP_TR153_TE58_etl61_0_8mm_gridding1_phc1_cfc0_sfc0_gnc0.json';
    %json_files{4} = 'D:\cartesian_maxgirf_epi_3d\data\long_phantom_20250222\meas_MID00035_FID10280_ep3d_tra_AP_TR111_TE37_etl37_0_8mm_trial1_gridding1_phc1_cfc0_sfc0_gnc0.json';
    %json_files{5} = 'D:\cartesian_maxgirf_epi_3d\data\long_phantom_20250222\meas_MID00046_FID10291_ep3d_tra_AP_highres_TR175_TE80_etl61_0_8mm_gridding1_phc1_cfc0_sfc0_gnc0.json';

    json_files{1} = 'D:\cartesian_maxgirf_epi_3d\data\long_phantom_20250222\meas_MID00035_FID10280_ep3d_tra_AP_TR111_TE37_etl37_0_8mm_trial1_gridding1_phc1_cfc0_sfc0_gnc0.json';
    %json_files{1} = 'D:\cartesian_maxgirf_epi_3d\data\long_phantom_20250222\meas_MID00035_FID10280_ep3d_tra_AP_TR111_TE37_etl37_0_8mm_trial2_gridding1_phc1_cfc0_sfc0_gnc0.json';
    %json_files{2} = 'D:\cartesian_maxgirf_epi_3d\data\long_phantom_20250222\meas_MID00035_FID10280_ep3d_tra_AP_TR111_TE37_etl37_0_8mm_trial3_gridding1_phc1_cfc0_sfc0_gnc0.json';
else
    json_files{1} = '/server/sdata/nlee/cartesian_maxgirf_epi_3d/data/long_phantom_20250222/meas_MID00025_FID10270_ep3d_tra_RL_TR217_TE90_etl101_gridding1_phc1_cfc0_sfc0_gnc0.json';
end

%% Calculate the number of json files
nr_json_files = length(json_files);  

%% Process per json file
for json_number = 1:nr_json_files

    %% Define the name of a .json file
    json_file = json_files{json_number};

    %% Calculate voxel coordinates
    demo_step1_cartesian_maxgirf_3d_calculate_voxel_coordinates;

    %% Calculate a fieldmap
    demo_calculate_fieldmap_gre_field_mapping;

    %% Prepare "imaging" k-space data
    demo_step3_cartesian_maxgirf_3d_prepare_ksp_imaging;

    %% Estimate CSMs
    demo_step4_cartesian_maxgirf_3d_estimate_csm;

    %% Cartesian MaxGIRF reconstruction
    demo_step5_cartesian_maxgirf_3d_recon;
end
