% demo_figure5.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/20/2024, Last modified: 03/09/2025

%% Clean slate
close all; clear all; clc;

%% Start a stopwatch timer
start_time = tic;

%% Set the full path of a directory
output_path1 = 'D:\cartesian_maxgirf_epi\data\GNL_MPRAGE_EPI\meas_MID00152_FID61451_ep2d_diff_4Trace_p2_fatsat_gridding1_phc1_conc1_gnl0_static0\PCLR';
output_path2 = 'D:\cartesian_maxgirf_epi\data\GNL_MPRAGE_EPI\meas_MID00152_FID61451_ep2d_diff_4Trace_p2_fatsat_gridding1_phc1_conc1_gnl1_static0\PCLR';

dicom_path1 = 'D:\cartesian_maxgirf_epi\data\GNL_MPRAGE_EPI\dicom\S25_ep2d_diff_4Trace_p2_fatsat_TRACEW';
dicom_path2 = 'D:\cartesian_maxgirf_epi\data\GNL_MPRAGE_EPI\dicom\S32_ep2d_diff_4Trace_p2_fatsat_TRACEW_S25_DIS2D';

%E:\projects_lenovo_20250319\cartesian_maxgirf_epi_2d\data

output_path1 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi\data\GNL_MPRAGE_EPI\meas_MID00152_FID61451_ep2d_diff_4Trace_p2_fatsat_gridding1_phc1_conc1_gnl0_static0\PCLR';
output_path2 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi\data\GNL_MPRAGE_EPI\meas_MID00152_FID61451_ep2d_diff_4Trace_p2_fatsat_gridding1_phc1_conc1_gnl1_static0\PCLR';

dicom_path1 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi\data\GNL_MPRAGE_EPI\dicom\S25_ep2d_diff_4Trace_p2_fatsat_TRACEW';
dicom_path2 = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi\data\GNL_MPRAGE_EPI\dicom\S32_ep2d_diff_4Trace_p2_fatsat_TRACEW_S25_DIS2D';

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

%--------------------------------------------------------------------------
% Number of slices
%--------------------------------------------------------------------------
N3 = nr_files / 2;

%% Read a .dicom file (Diffusion-weighted images without GNC)
img_dicom_gnc0 = zeros(N1, N2, N3, 'double');

for idx = N3+1:nr_files
    %----------------------------------------------------------------------
    % Get a DICOM header
    %----------------------------------------------------------------------
    dicom_info = dicominfo(fullfile(dir_info(idx).folder, dir_info(idx).name));

    %----------------------------------------------------------------------
    % Read a DICOM image
    %----------------------------------------------------------------------
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
    img_dicom_gnc0(:,:,idx-N3) = RescaleSlope * dicomread(dicom_info).' + RescaleIntercept; % transpose it!
end

%% Get directory information
dir_info = dir(fullfile(dicom_path2, '*IMA'));
nr_files = length(dir_info);

%% Read a .dicom file (Diffusion-weighted images with GNC)
img_dicom_gnc1 = zeros(N1, N2, N3, 'double');

for idx = N3+1:nr_files
    %----------------------------------------------------------------------
    % Get a DICOM header
    %----------------------------------------------------------------------
    dicom_info = dicominfo(fullfile(dir_info(idx).folder, dir_info(idx).name));

    %----------------------------------------------------------------------
    % Read a DICOM image
    %----------------------------------------------------------------------
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
    img_dicom_gnc1(:,:,idx-N3) = RescaleSlope * dicomread(dicom_info).' + RescaleIntercept; % transpose it!
end

%% Define variables
Nx = N1;
Ny = N2;

Nkx = 2 * Nx;
Nky = Ny;
nr_slices = nr_files / 2;
nr_diff_directions = 4;

%   "nr_slices": 20,
%   "nr_repetitions": 44,
%   "nr_diff_directions": 4,
%   "nr_diff_weightings": 2,
%   "nr_bvalues1": 4,
%   "nr_bvalues2": 10,
%   "main_orientation": 2

%% Read a .cfl file
img_maxgirf_gnc0 = complex(zeros(Nkx, Nky, nr_slices, nr_diff_directions, 'single')); % Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/0
img_maxgirf_gnc1 = complex(zeros(Nkx, Nky, nr_slices, nr_diff_directions, 'single')); % Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/1

for slice_number = 1:nr_slices
    %------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
    % slice_number: acquisition slice number
    % actual_slice_number is used to retrieve slice information from TWIX format
    %------------------------------------------------------------------
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

    for diff_dir_number = 1:nr_diff_directions
        %------------------------------------------------------------------
        % Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/0
        %------------------------------------------------------------------
        img1_filename1 = sprintf('img_maxgirf_PCLR_slc%d_wgt2_dir%d_gridding1_phc1_conc1_gnl0_static0_flat_i6_l0.00', slice_number, diff_dir_number);
        cfl_file = fullfile(output_path1, img1_filename1);
        tstart = tic; fprintf('%s:(SLC=%2d/%2d)(DIR=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices, diff_dir_number, nr_diff_directions);
        img_maxgirf_gnc0(:,:,actual_slice_number,diff_dir_number) = readcfl(cfl_file);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

        %------------------------------------------------------------------
        % Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/1
        %------------------------------------------------------------------
        img1_filename2 = sprintf('img_maxgirf_PCLR_slc%d_wgt2_dir%d_gridding1_phc1_conc1_gnl1_static0_flat_i6_l0.00', slice_number, diff_dir_number);
        cfl_file = fullfile(output_path2, img1_filename2);
        tstart = tic; fprintf('%s:(SLC=%2d/%2d)(DIR=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices, diff_dir_number, nr_diff_directions);
        img_maxgirf_gnc1(:,:,actual_slice_number,diff_dir_number) = readcfl(cfl_file);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end

%% Read a .cfl file
output_path = 'D:\cartesian_maxgirf_epi\data\GNL_MPRAGE_EPI\meas_MID00152_FID61451_ep2d_diff_4Trace_p2_fatsat_gridding1_phc1_conc1_gnl0_static0';

% Hack
output_path = 'E:\projects_lenovo_20250319\cartesian_maxgirf_epi\data\GNL_MPRAGE_EPI\meas_MID00152_FID61451_ep2d_diff_4Trace_p2_fatsat_gridding1_phc1_conc1_gnl0_static0';

slice_type = 'flat';

x = zeros(Nkx, Nky, nr_slices, 'single');
y = zeros(Nkx, Nky, nr_slices, 'single');
z = zeros(Nkx, Nky, nr_slices, 'single');

for slice_number = 1:nr_slices
    %------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
    % slice_number: acquisition slice number
    % actual_slice_number is used to retrieve slice information from TWIX format
    %------------------------------------------------------------------
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

    %----------------------------------------------------------------------
    % x (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('x_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    x(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('y_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    y(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('z_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    z(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Remove readout oversampling
idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(Nkx/2) + 1;

img_maxgirf_gnc0 = img_maxgirf_gnc0(idx1_range,:,:,:);
img_maxgirf_gnc1 = img_maxgirf_gnc1(idx1_range,:,:,:);

x = x(idx1_range,:,:);
y = y(idx1_range,:,:);
z = z(idx1_range,:,:);

%% Calculate trace images
img_maxgirf_gnc0_mean = mean(img_maxgirf_gnc0,4);
img_maxgirf_gnc1_mean = mean(img_maxgirf_gnc1,4);

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

FontSize = 14;
xlimits = [20 159];
ylimits = [20 159];

for slice_number = 9%1:nr_slices
    figure('Color', 'w', 'Position', [0 31 1128 947]);

    %----------------------------------------------------------------------
    % MaxGIRF + PCLR: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/0
    %----------------------------------------------------------------------
    ax1 = subplot(2,2,1);
    imagesc(abs(img_maxgirf_gnc1_mean(:,:,slice_number)).');
    colormap(ax1, gray(256));
    axis image off;
    xlim(ax1, xlimits);
    ylim(ax1, ylimits);
    text(ax1, N1 / 2, 12, '(Phase-const.) Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(ax1, N1 / 2, 19, {'corrections applied {\color[rgb]{0.9255 0.3882 0.0}during recon}'}, 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(ax1, 19, N2 / 2, {'w/ concomitant field correction', 'w/o GNL correction'}, 'Rotation', 90, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(ax1, N1 - 20, 4, '2D SE-EPI: axial, b1000, 1.45 x 1.45 mm^2, R = 2, PF = 6/8, ETL = 65, 10 NSA, mean DWI', 'Color', blue, 'Interpreter', 'tex', 'FontSize', FontSize + 2, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(ax1, 21, 19, '(A)', 'FontSize', 16, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(ax1, N2 / 2, ylimits(1), 'complex averaging', 'FontSize', FontSize, 'Rotation', 0, 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

    text(ax1, ylimits(1) - 1, xlimits(2) + 5, {'Slice at', sprintf('z = %4.2f mm', z(1,1,slice_number) * 1e3)}, 'FontSize', FontSize, 'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    text(ax1, ylimits(1), N2 / 2, sprintf('PE direction (A >> P)'), 'FontSize', FontSize, 'Rotation', -90, 'Interpreter', 'tex', 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
        
    %----------------------------------------------------------------------
    % MaxGIRF + PCLR: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/1
    %----------------------------------------------------------------------
    ax3 = subplot(2,2,3);
    imagesc(abs(img_maxgirf_gnc0_mean(:,:,slice_number)).');
    colormap(ax3, gray(256));
    axis image off;
    xlim(ax3, xlimits);
    ylim(ax3, ylimits);
    text(ax3, 19, N2 / 2, {'w/ concomitant field correction', 'w/ GNL correction'}, 'Rotation', 90, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(ax3, 21, 19, '(C)', 'FontSize', 16, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(ax3, N2 / 2, ylimits(1), 'complex averaging', 'FontSize', FontSize, 'Rotation', 0, 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

    %----------------------------------------------------------------------
    % DICOM: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/0
    %----------------------------------------------------------------------
    ax2 = subplot(2,2,2);
    imagesc(img_dicom_gnc0(:,:,slice_number).');
    colormap(ax2, gray(256));
    axis image off;
    xlim(ax2, xlimits);
    ylim(ax2, ylimits);
    text(ax2, N1 / 2, 12, 'Traditional reconstruction', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(ax2, N1 / 2, 19, 'corrections applied {\color[rgb]{0.0118 0.6000 0.6000}after recon}', 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(ax2, 21, 19, '(B)', 'FontSize', 16, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(ax2, N2 / 2, ylimits(1), 'magnitude averaging', 'FontSize', FontSize, 'Rotation', 0, 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

    %----------------------------------------------------------------------
    % DICOM: Gridding/PHC/CFC/SFC/GNC = 1/1/1/0/1
    %----------------------------------------------------------------------
    ax4 = subplot(2,2,4);
    imagesc(img_dicom_gnc1(:,:,slice_number).');
    colormap(ax4, gray(256));
    axis image off;
    xlim(ax4, xlimits);
    ylim(ax4, ylimits);
    text(ax4, 21, 19, '(D)', 'FontSize', 16, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    text(ax4, N2 / 2, ylimits(1), 'magnitude averaging', 'FontSize', FontSize, 'Rotation', 0, 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

    set(ax1, 'Position', [0.1516 0.4613 0.4124 0.4468]);
    set(ax2, 'Position', [0.5280 0.4613 0.4124 0.4468]);
    set(ax3, 'Position', [0.1516 0.0129 0.4124 0.4468]);
    set(ax4, 'Position', [0.5280 0.0129 0.4124 0.4468]);

    %----------------------------------------------------------------------
    % Arrow 1 (left)
    %----------------------------------------------------------------------
    % Create arrow

    d1 = 0.220744680851063 - 0.245567375886524;
    d2 = 0.240248226950354 - 0.265070921985815;
    d3 = 0.575501583949309 - 0.591341077085532;
    d4 = 0.590285110876443 - 0.606124604012666;

    annotation(gcf, 'arrow', [0.245567375886524+d1 0.265070921985815+d2],...
        [0.591341077085532+d3 0.606124604012666+d4], 'Color', [1 0 0], 'LineWidth', 2,...
        'HeadStyle', 'plain');
    
    % Create arrow
    annotation(gcf, 'arrow', [0.62234042553191+d1 0.641843971631201+d2],...
        [0.591341077085532+d3 0.606124604012666+d4], 'Color', [1 0 0], 'LineWidth', 2,...
        'HeadStyle', 'plain');

    % Create arrow
    annotation(gcf, 'arrow', [0.220744680851063 0.240248226950354],...
        [0.12671594508974 0.141499472016874], 'Color', [1 0 0], 'LineWidth', 2,...
        'HeadStyle', 'plain');

    % Create arrow
    annotation(gcf, 'arrow', [0.597517730496449 0.61702127659574],...
        [0.129883843716984 0.144667370644118], 'Color', [1 0 0], 'LineWidth', 2,...
        'HeadStyle', 'plain');

    %----------------------------------------------------------------------
    % Arrow 1 (right)
    %----------------------------------------------------------------------
    % Create arrow
    annotation(gcf, 'arrow', [0.398936170212766 0.382092198581558],...
        [0.682154171066526 0.665258711721214], 'Color', color_order{3}, 'LineWidth', 2,...
        'HeadStyle', 'plain');

    % Create arrow
    annotation(gcf, 'arrow', [0.775709219858151 0.758865248226944],...
        [0.682154171066526 0.665258711721214], 'Color', color_order{3}, 'LineWidth', 2,...
        'HeadStyle', 'plain');

    % Create arrow
    annotation(gcf, 'arrow', [0.398936170212766 0.382092198581558],...
        [0.229144667370628 0.212249208025316], 'Color', color_order{3}, 'LineWidth', 2,...
        'HeadStyle', 'plain');

    % Create arrow
    annotation(gcf, 'arrow', [0.775709219858151 0.758865248226944],...
        [0.229144667370628 0.212249208025316], 'Color', color_order{3}, 'LineWidth', 2,...
        'HeadStyle', 'plain');

    export_fig(sprintf('figure5_slc%d', slice_number), '-r300', '-tif', '-c[20, 260, 24, 450]'); % [top,right,bottom,left]
    close gcf;
end
