% demo_cartesian_maxgirf_2d_recon_lowres.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/18/2025, Last modified: 02/08/2025

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
Lmax          = json.recon_parameters.Lmax;          % maximum rank of the SVD approximation of a higher-order encoding matrix
L             = json.recon_parameters.L;             % rank of the SVD approximation of a higher-order encoding matrix
lambda        = json.recon_parameters.lambda;        % l2 regularization parameter
tol           = json.recon_parameters.tol;           % PCG tolerance
maxiter       = json.recon_parameters.maxiter;       % PCG maximum iteration 
slice_type    = json.recon_parameters.slice_type;    % type of an excitation slice: "curved" vs "flat"
cal_size      = json.recon_parameters.cal_size.';    % size of calibration region
phc_flag      = json.recon_parameters.phc_flag;      % 1=yes, 0=no
gridding_flag = json.recon_parameters.gridding_flag; % 1=yes, 0=no
cfc_flag      = json.recon_parameters.cfc_flag;      % 1=yes, 0=no
sfc_flag      = json.recon_parameters.sfc_flag;      % 1=yes, 0=no
gnc_flag      = json.recon_parameters.gnc_flag;      % 1=yes, 0=no

%--------------------------------------------------------------------------
% Number of slices
%--------------------------------------------------------------------------
if isfield(json, 'nr_slices')
    nr_slices = json.nr_slices;
else
    nr_slices = 1;
end

%--------------------------------------------------------------------------
% Number of repetitions
%--------------------------------------------------------------------------
if isfield(json, 'nr_repetitions')
    nr_repetitions = json.nr_repetitions;
else
    nr_repetitions = 1;
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

%% Calculate the number of total slices
nr_recons = nr_slices * nr_repetitions;

%% Perform image reconstruction per slice
for idx = 1:nr_recons

    %% Get information about the current slice
    [slice_number, repetition_number] = ind2sub([nr_slices nr_repetitions], idx);

    % 42
    if slice_number ~= 13
        continue;
    end

    %% Read a .cfl file
    %----------------------------------------------------------------------
    % ksp_cartesian (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('ksp_img_cartesian_slc%d_rep%d_gridding%d_phc%d', slice_number, repetition_number, gridding_flag, phc_flag));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    ksp = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % mask_img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('mask_img_slc%d', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    img_mask = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % mask_cal (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('mask_cal_slc%d', slice_number));
    if exist(strcat(cfl_file, '.cfl'), 'file') % R > 1
        tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
        cal_mask = readcfl(cfl_file);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    else % either R=1 or PF
        cal_mask = img_mask .* flip(img_mask,2);
    end

    %----------------------------------------------------------------------
    % circle_mask (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('circle_mask_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    circle_mask = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % sens (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('sens_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    sens = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % ev_maps (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('ev_maps_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    ev_maps = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % x (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('x_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    x = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('y_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    y = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('z_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    z = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % u (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('u_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    u = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % v (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('v_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    v = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % w (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('w_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    w = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % du (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('du_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    du = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dv (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dv_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    dv = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dw (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dw_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    dw = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % read_sign
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, 'read_sign');
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    read_sign = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % phase_sign
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, 'phase_sign');
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    phase_sign = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % initial_phase (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('initial_phase_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    initial_phase = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % U (NkNs x Lmax)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('U_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    U = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % S (Lmax x 1)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('S_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    S = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % V (Nkx x Nky x Nkz x Lmax)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('V_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    V = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Apply a threshold mask on CSMs
%     % 0.98, 5 => 0.99, 10
%     ev_mask = (ev_maps > 0.96); % 0.94 for in-vivo
%     ev_mask2 = bwareaopen(ev_mask, 600); % Keep only blobs with an area of 300 pixels or more.
%     se = strel('disk',5);
%     ev_mask_dilated = imdilate(ev_mask2,se);
%     ev_mask_dilated = imfill(ev_mask_dilated, "holes");
%     sens = bsxfun(@times, sens, ev_mask_dilated);

    ev_mask = (ev_maps > 0.99); % 0.94 for in-vivo
    se = strel('disk',20);
    ev_mask_opened = imopen(ev_mask,se);
    se = strel('disk',5);
    ev_mask_dilated = imdilate(ev_mask_opened,se);
    sens = bsxfun(@times, sens, ev_mask_dilated);

    %%
    if 0
    figure;
    subplot(4,1,1);
    imagesc(ev_mask.');
    axis image;

    subplot(4,1,2);
    imagesc(ev_mask_opened.');
    axis image;

    subplot(4,1,3);
    imagesc(ev_mask_dilated.');
    axis image;
    end
    
    
    %% Calculate a low-resolution sampling mask
    lowres_mask = img_mask .* cal_mask;

    %% Subselect U
    img_list = find(img_mask);
    lowres_list = find(lowres_mask);
    U_lowres = U(ismember(img_list, lowres_list),:);

    %% Calculate low-resolution only k-space data
    ksp_lowres = bsxfun(@times, lowres_mask, ksp);

    %% Calculate spatial positions for Type-1 and Type-2 NUFFTs
    if gnc_flag
        p1 = (u + du) / encoded_fov(1) * (2 * pi); % RO [-pi,pi]
        p2 = (v + dv) / encoded_fov(2) * (2 * pi); % PE [-pi,pi]
    else
        p1 = u / encoded_fov(1) * (2 * pi); % RO [-pi,pi]
        p2 = v / encoded_fov(2) * (2 * pi); % PE [-pi,pi]
    end

    %% Calculate a support mask
    support_mask = zeros(Nkx, Nky, Nkz, 'single');
    support_mask((abs(p1) < pi) & (abs(p2) < pi) & (circle_mask > 0)) = 1;

    %% Select spatial positions within a support mask
    p1 = p1(support_mask > 0);
    p2 = p2(support_mask > 0);

    %% Set parameters for Type-1 and Type-2 NUFFTs
    eps = 1e-6;
    iflag = -1;

    %% Type-1 NUFFT based Cartesian MaxGIRF operators
    Ah = @(x) cartesian_maxgirf_2d_adjoint(x, sens, lowres_mask, p1, p2, iflag, eps, support_mask, U_lowres, V, L, initial_phase);
    AhA = @(x) cartesian_maxgirf_2d_normal(x, sens, lowres_mask, p1, p2, iflag, eps, support_mask, U_lowres, V, L, initial_phase, lambda);

    %% Perform Cartesian MaxGIRF reconstruction
    clear cartesian_maxgirf_2d_normal;
    tstart = tic; fprintf('%s:(SLC=%2d/%2d)(REP=%2d/%2d) Performing Cartesian MaxGIRF reconstruction (low-resolution):\n', datetime, slice_number, nr_slices, repetition_number, nr_repetitions);
    [img, flag, relres, iter, resvec] = pcg(@(x) AhA(x), Ah(ksp_lowres), tol, maxiter);
    img = reshape(img, [Nkx Nky Nkz]);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('img_maxgirf_lowres_slc%d_rep%d_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_%s_i%d_l%4.2f', slice_number, repetition_number, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, slice_type, maxiter, lambda);
    cfl_file = fullfile(output_path, img_filename);
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, img);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % support_mask (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    mask_filename = sprintf('support_mask_slc%d_%s', slice_number, slice_type);
    cfl_file = fullfile(output_path, mask_filename);
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, support_mask);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Display the magnitude of an image
    FontSize = 14;

    c1 = floor(Nkx/2) + 1;
    c2 = floor(Nky/2) + 1;
    c3 = floor(Nkz/2) + 1;

    xmax = max([abs(max(x(:))) abs(min(x(:)))]);
    ymax = max([abs(max(y(:))) abs(min(y(:)))]);
    zmax = max([abs(max(z(:))) abs(min(z(:)))]);

    xlimits = [-xmax xmax];
    ylimits = [-ymax ymax];
    zlimits = [-zmax zmax];

    if main_orientation == 2 % TRANSVERSAL = 2
        Position = [1020 44 525 934];
        slice_direction = 'z';
        slice_offset = z(c1,c2) * 1e3;
        ax = 0;
        el = 90;
    elseif main_orientation == 0 % SAGITTAL = 0
        Position = [1038 25 778 953];
        slice_direction = 'x';
        slice_offset = x(c1,c2) * 1e3;
        ax = -90;
        el = 0;
    elseif main_orientation == 1 % CORONAL = 1
        Position = [680 28 1029 950];
        slice_direction = 'y';
        slice_offset = y(c1,c2) * 1e3;
        ax = 0;
        el = 0;
    end

    if read_sign == -1
        x = flip(x,1);
        y = flip(y,1);
        z = flip(z,1);
        img = flip(img,1);
    end

    if phase_sign == -1
        x = flip(x,2);
        y = flip(y,2);
        z = flip(z,2);
        img = flip(img,2);
    end

    figure('Color', 'w', 'Position', Position);
    surf(x * 1e3, y * 1e3, z * 1e3, abs(img), 'EdgeColor', 'none');
    colormap(gray(256));
    axis image;
    xlim(xlimits * 1e3);
    ylim(ylimits * 1e3);
    zlim(zlimits * 1e3);
    %caxis([0 15]);

    title_text1 = sprintf('SLC/REP = %d/%d, %s slice, %s = %4.1f mm', slice_number, repetition_number, slice_type, slice_direction, slice_offset);
    title_text2 = sprintf('Gridding/PHC/CFC/SFC/GNC = %d/%d/%d/%d/%d', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag);
    title_text3 = sprintf('max. iterations = %d, $$\\lambda$$ = %4.2f', maxiter, lambda);
    fig_filename = sprintf('img_maxgirf_lowres_slc%d_rep%d_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_%s_i%d_l%4.2f', slice_number, repetition_number, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, slice_type, maxiter, lambda);

    title({'Cartesian MaxGIRF (low-resolution)', title_text1, title_text2, title_text3}, 'FontWeight', 'normal', 'FontSize', FontSize, 'Interpreter', 'latex');
    xlabel('x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel('y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel('z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    set(gca, 'TickLabelInterpreter', 'latex');
    view(ax,el);
    set(gca, 'ZDir', 'reverse');
    drawnow;
    export_fig(fullfile(output_path, sprintf('%s_mag', fig_filename)), '-r300', '-tif');
    close gcf;

    %% Display the phase of an image
    figure('Color', 'w', 'Position', Position);
    surf(x * 1e3, y * 1e3, z * 1e3, angle(img) * 180 / pi, 'EdgeColor', 'none');
    colormap(hsv(256));
    axis image;
    xlim(xlimits * 1e3);
    ylim(ylimits * 1e3);
    zlim(zlimits * 1e3);
    caxis([-180 180]);
    title({'Cartesian MaxGIRF (low-resolution)', title_text1, title_text2, title_text3}, 'FontWeight', 'normal', 'FontSize', FontSize, 'Interpreter', 'latex');
    xlabel('x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel('y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel('z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    set(gca, 'TickLabelInterpreter', 'latex');
    view(ax,el);
    set(gca, 'ZDir', 'reverse');
    drawnow;
    export_fig(fullfile(output_path, sprintf('%s_phase', fig_filename)), '-r300', '-tif');
    close gcf;

    %% Display a support mask
    figure('Color', 'w', 'Position', Position);
    surf(x * 1e3, y * 1e3, z * 1e3, support_mask, 'EdgeColor', 'none');
    colormap(gray(256));
    axis image;
    xlim(xlimits * 1e3);
    ylim(ylimits * 1e3);
    zlim(zlimits * 1e3);
    caxis([0 1]);
    title({'Cartesian MaxGIRF (mask)', sprintf('SLC = %d, %s slice, %s = %4.1f mm', slice_number, slice_type, slice_direction, slice_offset), ...
        title_text2, title_text3}, ...
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
end
