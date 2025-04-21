% batch_phantom0331_20240828_FH_vs_HF.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 11/26/2024, Last modified: 01/31/2025

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
    %json_files{1} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00063_FID21492_ep2d_se_bw1002_cor_FH_gridding1_phc1_cfc0_sfc0_gnc0.json';
    %json_files{2} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00063_FID21492_ep2d_se_bw1002_cor_FH_gridding1_phc1_cfc1_sfc0_gnc0.json';
    %json_files{3} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00063_FID21492_ep2d_se_bw1002_cor_FH_gridding1_phc1_cfc1_sfc0_gnc1.json';
    json_files{1} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00063_FID21492_ep2d_se_bw1002_cor_FH_gridding1_phc1_cfc1_sfc1_gnc1.json';

    %json_files{1} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00065_FID21494_ep2d_se_bw1002_cor_HF_gridding1_phc1_cfc0_sfc0_gnc0.json';
    %json_files{2} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00065_FID21494_ep2d_se_bw1002_cor_HF_gridding1_phc1_cfc1_sfc0_gnc0.json';
    %json_files{3} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00065_FID21494_ep2d_se_bw1002_cor_HF_gridding1_phc1_cfc1_sfc0_gnc1.json';
    %json_files{1} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00065_FID21494_ep2d_se_bw1002_cor_HF_gridding1_phc1_cfc1_sfc1_gnc1.json';
else
    json_files{1} = '/server/sdata/nlee/cartesian_maxgirf_epi_2d/data/phantom0331_20240828/meas_MID00063_FID21492_ep2d_se_bw1002_cor_FH_gridding1_phc1_cfc0_sfc0_gnc0.json';
end

%% Define input files for staic field correction (static off-resonance correction)
flash_json_files{1} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00068_FID21497_flash_te37_gnc1.json';
flash_json_files{2} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00069_FID21498_flash_te47_gnc1.json';
flash_json_files{3} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00070_FID21499_flash_te57_gnc1.json';
flash_json_files{4} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00071_FID21500_flash_te67_gnc1.json';
flash_json_files{5} = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00072_FID21501_flash_te77_gnc1.json';

%% Calculate the number of json files
nr_json_files = length(json_files);

%% Process per json file
for json_number = 1:nr_json_files

    %% Define the name of a .json file
    json_file = json_files{json_number};

    %% Calculate voxel coordinates
    demo_step1_cartesian_maxgirf_2d_calculate_voxel_coordinates;

    %% Calculate a fieldmap
    demo_calculate_fieldmap_single_echo_flash_2d;

    %% Prepare "calibration" k-space data
    %demo_step2_cartesian_maxgirf_prepare_ksp_calibration;

    %% Prepare "imaging" k-space data
    demo_step3_cartesian_maxgirf_2d_prepare_ksp_imaging;

    %% Estimate CSMs
    demo_step4_cartesian_maxgirf_2d_estimate_csm;

    %% Cartesian MaxGIRF reconstruction
    demo_step5_cartesian_maxgirf_2d_recon;

    %% Cartesian MaxGIRF reconstruction (low-resolution)
    %demo_step6_cartesian_maxgirf_recon_low_resolution;

    %% Phase-constrained Cartesian MaxGIRF & low-rank and subspace model-based reconstruction
    %demo_step7_cartesian_maxgirf_recon_phase_constrained_lowrank;
end
