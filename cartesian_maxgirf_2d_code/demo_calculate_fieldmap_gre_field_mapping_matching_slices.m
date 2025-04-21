% demo_calculate_fieldmap_gre_field_mapping_matching_slices.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/20/2025, Last modified: 02/12/2025

%% Clean slate
close all; clearvars -except json_number nr_json_files json_files json_file grad_file_path fieldmap_json_file; clc;

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
    output_path       = strrep(json.output_path, '/', '\');
    ismrmrd_data_file = strrep(json.ismrmrd_data_file, '/', '\');
else
    output_path       = json.output_path;
    ismrmrd_data_file = json.ismrmrd_data_file;
end

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
slice_type = json.recon_parameters.slice_type; % type of an excitation slice: "curved" vs "flat"
sfc_flag   = json.recon_parameters.sfc_flag;   % 1=yes, 0=no

%--------------------------------------------------------------------------
% Number of slices
%--------------------------------------------------------------------------
if isfield(json, 'nr_slices')
    nr_slices = json.nr_slices;
else
    nr_slices = 1;
end

%% Check if static field correction is on
if sfc_flag == 0
    return;
end

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
Nkx = header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky = header.encoding.encodedSpace.matrixSize.y; % number of phase-encoding steps in k-space
Nkz = header.encoding.encodedSpace.matrixSize.z; % number of slice-encoding steps in k-space

%--------------------------------------------------------------------------
% Recon Space (Nx, Ny, Nz)
%--------------------------------------------------------------------------
Nx = header.encoding.reconSpace.matrixSize.x; % number of samples in image space (RO)
Ny = header.encoding.reconSpace.matrixSize.y; % number of samples in image space (PE)
Nz = header.encoding.reconSpace.matrixSize.z; % number of samples in image space (SL)

%% Load a .cfi file
x = zeros(Nkx, Nky, nr_slices, 'single');
y = zeros(Nkx, Nky, nr_slices, 'single');
z = zeros(Nkx, Nky, nr_slices, 'single');

for slice_number = 1:nr_slices
    %----------------------------------------------------------------------
    % x (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('x_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    x(:,:,slice_number) = real(readcfl(cfl_file));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('y_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    y(:,:,slice_number) = real(readcfl(cfl_file));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('z_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    z(:,:,slice_number) = real(readcfl(cfl_file));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, fieldmap_json_file);
fid_fieldmap = fopen(fieldmap_json_file); 
json_txt_fieldmap = fread(fid_fieldmap, [1 inf], 'char=>char'); 
fclose(fid_fieldmap);
json_fieldmap = jsondecode(json_txt_fieldmap);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
if ispc
    siemens_twix_file_fieldmap = strrep(json.siemens_twix_file, '/', '\');
    output_path_fieldmap       = strrep(json_fieldmap.output_path, '/', '\');
    ismrmrd_data_file_fieldmap = strrep(json_fieldmap.ismrmrd_data_file, '/', '\');
else
    siemens_twix_file_fieldmap = json.siemens_twix_file;
    output_path_fieldmap       = json_fieldmap.output_path;
    ismrmrd_data_file_fieldmap = json_fieldmap.ismrmrd_data_file;
end

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
lambda     = json_fieldmap.recon_parameters.lambda;     % l2 regularization parameter
tol        = json_fieldmap.recon_parameters.tol;        % PCG tolerance
maxiter    = json_fieldmap.recon_parameters.maxiter;    % PCG maximum iteration 
slice_type = json_fieldmap.recon_parameters.slice_type; % type of an excitation slice: "curved" vs "flat"
gnc_flag   = json_fieldmap.recon_parameters.gnc_flag;   % 1=yes, 0=no

%--------------------------------------------------------------------------
% Number of slices
%--------------------------------------------------------------------------
if isfield(json_fieldmap, 'nr_slices')
    nr_slices_fieldmap = json_fieldmap.nr_slices;
else
    nr_slices_fieldmap = 1;
end

%--------------------------------------------------------------------------
% Smoothing parameters
%--------------------------------------------------------------------------
if isfield(json_fieldmap, 'smoothing_parameters')
    w_fB0 = json_fieldmap.smoothing_parameters.w_fB0; % 32 in BART?
    h_fB0 = json_fieldmap.smoothing_parameters.h_fB0; % Sobolev index for B0 field inhomogeneity
else
    w_fB0 = 32;
    h_fB0 = 2;
end

h_fB0 = 4;

%% Read an ISMRMRD file (k-space data)
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... ', datetime, ismrmrd_data_file_fieldmap);
if exist(ismrmrd_data_file_fieldmap, 'file')
    dset = ismrmrd.Dataset(ismrmrd_data_file_fieldmap, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_file_fieldmap);
end

%% Get imaging parameters from an XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% measurement information
%--------------------------------------------------------------------------
patient_position = header.measurementInformation.patientPosition;

%--------------------------------------------------------------------------
% Encoding Space (Nkx, Nky, Nkz)
%--------------------------------------------------------------------------
encoded_fov_fieldmap(1) = header.encoding.encodedSpace.fieldOfView_mm.x * 1e-3; % [m] RO
encoded_fov_fieldmap(2) = header.encoding.encodedSpace.fieldOfView_mm.y * 1e-3; % [m] PE
encoded_fov_fieldmap(3) = header.encoding.encodedSpace.fieldOfView_mm.z * 1e-3; % [m] SL

Nkx_fieldmap = header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky_fieldmap = header.encoding.encodedSpace.matrixSize.y; % number of phase-encoding steps in k-space
Nkz_fieldmap = header.encoding.encodedSpace.matrixSize.z; % number of slice-encoding steps in k-space

encoded_resolution_fieldmap = encoded_fov_fieldmap ./ [Nkx_fieldmap Nky_fieldmap Nkz_fieldmap]; % [m]

%--------------------------------------------------------------------------
% Recon Space (Nx, Ny, Nz)
%--------------------------------------------------------------------------
Nx_fieldmap = header.encoding.reconSpace.matrixSize.x; % number of samples in image space (RO)
Ny_fieldmap = header.encoding.reconSpace.matrixSize.y; % number of samples in image space (PE)
Nz_fieldmap = header.encoding.reconSpace.matrixSize.z; % number of samples in image space (SL)

%% Read acquisitions
tstart = tic; fprintf('%s: Reading acquisitions... ', datetime);
raw_data = dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Get data type
%--------------------------------------------------------------------------
acq_is_noise_measurement                = raw_data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
acq_is_parallel_calibration             = raw_data.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION');
acq_is_parallel_calibration_and_imaging = raw_data.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING');
acq_is_reverse                          = raw_data.head.flagIsSet('ACQ_IS_REVERSE');
acq_is_navigation_data                  = raw_data.head.flagIsSet('ACQ_IS_NAVIGATION_DATA');
acq_is_phasecorr_data                   = raw_data.head.flagIsSet('ACQ_IS_PHASECORR_DATA');
acq_is_hpfeedback_data                  = raw_data.head.flagIsSet('ACQ_IS_HPFEEDBACK_DATA');
acq_is_dummyscan_data                   = raw_data.head.flagIsSet('ACQ_IS_DUMMYSCAN_DATA');
acq_is_rtfeedback_data                  = raw_data.head.flagIsSet('ACQ_IS_RTFEEDBACK_DATA');
acq_is_surfacecoilcorrectionscan_data   = raw_data.head.flagIsSet('ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA');

%% Get imaging data
img_data = raw_data.select(find(~acq_is_noise_measurement & ~acq_is_parallel_calibration & ~acq_is_navigation_data & ~acq_is_phasecorr_data));

%% Load a .cfl file
%--------------------------------------------------------------------------
% TE (N1 x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path_fieldmap, 'TE');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
TE = real(readcfl(cfl_file)).';
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Load a .cfl file
img_fieldmap = zeros([Nkx_fieldmap Nky_fieldmap nr_slices_fieldmap 1 1 2], 'single');  % N1 x N2 x N3 x 1 x 1 x Ne

x_fieldmap = zeros(Nkx_fieldmap, Nky_fieldmap, nr_slices_fieldmap, 'single');
y_fieldmap = zeros(Nkx_fieldmap, Nky_fieldmap, nr_slices_fieldmap, 'single');
z_fieldmap = zeros(Nkx_fieldmap, Nky_fieldmap, nr_slices_fieldmap, 'single');

slice_number_list = zeros(nr_slices_fieldmap, 1, 'single');

for slice_number = 1:nr_slices_fieldmap
    %----------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
    %----------------------------------------------------------------------
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
    slice_number_list(actual_slice_number) = slice_number;

    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('img_type1_slc%d_eco1_gnc%d_%s_i%d_l%4.2f', slice_number, gnc_flag, slice_type, maxiter, lambda);
    cfl_file = fullfile(output_path_fieldmap, img_filename);
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    img_fieldmap(:,:,actual_slice_number,1,1,1) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('img_type1_slc%d_eco2_gnc%d_%s_i%d_l%4.2f', slice_number, gnc_flag, slice_type, maxiter, lambda);
    cfl_file = fullfile(output_path_fieldmap, img_filename);
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    img_fieldmap(:,:,actual_slice_number,1,1,2) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % x (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path_fieldmap, sprintf('x_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    x_fieldmap(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path_fieldmap, sprintf('y_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    y_fieldmap(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path_fieldmap, sprintf('z_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    z_fieldmap(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%--------------------------------------------------------------------------
% read_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'read_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
read_sign_fieldmap = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% phase_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'phase_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
phase_sign_fieldmap = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a Siemens .dat file
fprintf('%s: Reading a Siemens .dat file: %s\n', datetime, siemens_twix_file_fieldmap);
twix = mapVBVD(siemens_twix_file_fieldmap);

%% Get the name of a gradient set
if isfield(twix{1}.hdr.Meas, 'tGradientCoil')
    tGradientCoil = twix{1}.hdr.Meas.tGradientCoil;
elseif isfield(twix{2}.hdr.Meas, 'tGradientCoil')
    tGradientCoil = twix{2}.hdr.Meas.tGradientCoil;
end

%% Reduce the TWIX dataset
if length(twix) > 1
    twix = twix{end};
end

%% Get a rotation matrix from the GCS to the PCS (ISMRMRD format)
phase_dir = double(img_data.head.phase_dir(:,1));
read_dir  = double(img_data.head.read_dir(:,1));
slice_dir = double(img_data.head.slice_dir(:,1));
R_gcs2pcs_ismrmrd = [phase_dir read_dir slice_dir];

%% Calculate a transformation matrix from the RCS to the GCS [r,c,s] <=> [PE,RO,SL]
R_rcs2gcs = [0    1    0 ; % [PE]   [0 1 0] * [r]
             1    0    0 ; % [RO] = [1 0 0] * [c]
             0    0    1]; % [SL]   [0 0 1] * [s]

%% Calculate a rotation matrix from the PCS to the DCS
R_pcs2dcs = siemens_calculate_transform_pcs_to_dcs(patient_position);

%% Calculate a rotation matrix from the GCS to the DCS
R_gcs2dcs = R_pcs2dcs * R_gcs2pcs_ismrmrd;

%% Calculate a scaling matrix [m]
scaling_matrix = diag(encoded_resolution_fieldmap / 2);

%% Calculate spatial coordinates in the RCS [m]
if read_sign_fieldmap < 0
    row_range = ((1:Nkx).' - (floor(Nkx/2) + 1)) + 1;
else
    row_range = ((1:Nkx).' - (floor(Nkx/2) + 1));
end
 
if phase_sign_fieldmap < 0
    col_range = ((1:Nky).' - (floor(Nky/2) + 1)) + 1;
else
    col_range = ((1:Nky).' - (floor(Nky/2) + 1));
end

slc_range = (1:Nkz).' - (floor(Nkz/2) + 1);

[I1,I2,I3] = ndgrid(row_range, col_range, slc_range);
r_rcs = scaling_matrix * cat(1, I1(:).', I2(:).', I3(:).'); % 3 x N

%% Calculate ideal spatial coordinates in the GCS [m]
r_gcs = R_rcs2gcs * r_rcs; % 3 x N

%% Calculate voxel coordinates
x_fieldmap_interp = zeros(2 * Nkx_fieldmap, 2 * Nky_fieldmap, nr_slices_fieldmap, 'single');
y_fieldmap_interp = zeros(2 * Nkx_fieldmap, 2 * Nky_fieldmap, nr_slices_fieldmap, 'single');
z_fieldmap_interp = zeros(2 * Nkx_fieldmap, 2 * Nky_fieldmap, nr_slices_fieldmap, 'single');

for idx = 1:nr_slices

    %% Get information about the current slice
    slice_number = ind2sub(nr_slices, idx);

    %% Get a list of imaging acquisitions
    img_acq_list = find((img_data.head.idx.slice == (slice_number - 1)) & (img_data.head.idx.repetition == 0));

    %% Calculate the actual slice number for Siemens interleaved multislice imaging
    %----------------------------------------------------------------------
    % slice_number: acquisition slice number
    % actual_slice_number is used to retrieve slice information from TWIX format
    % Q: For coronal? acq_slice_nr = slice_nr;
    %----------------------------------------------------------------------
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

    %% Get a slice offset in the PCS from Siemens TWIX format
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}, 'sPosition')
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition, 'dSag')
            sag_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition.dSag; % [mm]
        else
            sag_offset_twix = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition, 'dCor')
            cor_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition.dCor; % [mm]
        else
            cor_offset_twix = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition, 'dTra')
            tra_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition.dTra; % [mm]
        else
            tra_offset_twix = 0; % [mm]
        end
    else
        sag_offset_twix = 0; % [mm]
        cor_offset_twix = 0; % [mm]
        tra_offset_twix = 0; % [mm]
    end

    %% Get a slice offset of a stack in the PCS from Siemens TWIX format
    pcs_offset = [sag_offset_twix; cor_offset_twix; tra_offset_twix] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Get a slice offset in the PCS from ISMRMRD format
    sag_offset_ismrmrd = double(img_data.head.position(1,img_acq_list(1))); % [mm]
    cor_offset_ismrmrd = double(img_data.head.position(2,img_acq_list(1))); % [mm]
    tra_offset_ismrmrd = double(img_data.head.position(3,img_acq_list(1))); % [mm]
    pcs_offset_ismrmrd = [sag_offset_ismrmrd; cor_offset_ismrmrd; tra_offset_ismrmrd] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Calculate a slice offset in the DCS [m]
    dcs_offset = R_pcs2dcs * pcs_offset; % 3 x 1

    %% Calculate ideal spatial coordinates in the DCS [m]
    r_dcs = R_gcs2dcs * r_gcs + repmat(dcs_offset, [1 2 * Nkx_fieldmap * 2 * Nky_fieldmap * Nkz_fieldmap]); % 3 x N
    x_fieldmap_interp(:,:,actual_slice_number) = reshape(r_dcs(1,:), [2 * Nkx_fieldmap 2 * Nky_fieldmap Nkz_fieldmap]); % x [m]
    y_fieldmap_interp(:,:,actual_slice_number) = reshape(r_dcs(2,:), [2 * Nkx_fieldmap 2 * Nky_fieldmap Nkz_fieldmap]); % y [m]
    z_fieldmap_interp(:,:,actual_slice_number) = reshape(r_dcs(3,:), [2 * Nkx_fieldmap 2 * Nky_fieldmap Nkz_fieldmap]); % z [m]
end

%% Zeropad in k-space (fieldmap)
idx1_range = (-floor(Nkx_fieldmap/2):ceil(Nkx_fieldmap/2)-1).' + floor(2 * Nkx_fieldmap / 2) + 1;
idx2_range = (-floor(Nky_fieldmap/2):ceil(Nky_fieldmap/2)-1).' + floor(2 * Nky_fieldmap / 2) + 1;

ksp_fieldmap_zpad = complex(zeros([2 * Nkx_fieldmap 2 * Nky_fieldmap nr_slices 1 1 2], 'single'));
for echo_number = 1:2
    %--------------------------------------------------------------------------
    % Siemens: k-space <=> image space
    %--------------------------------------------------------------------------
    tstart = tic; fprintf('%s: Applying inverse FFT to move from image space to k-space... ', datetime);
    ksp_fieldmap = img_fieldmap(:,:,:,1,1,echo_number);
    for dim = 1:2
        ksp_fieldmap = sqrt(size(ksp_fieldmap,dim)) * fftshift(ifft(ifftshift(ksp_fieldmap, dim), [], dim), dim);
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    ksp_fieldmap_zpad(idx1_range, idx2_range, :, :, :, echo_number) = ksp_fieldmap;
end

%% Perform sinc interpolation in image space
%--------------------------------------------------------------------------
% Siemens: k-space <=> image space
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Applying forward FFT to move from k-space to image space... ', datetime);
img_fieldmap_interp = ksp_fieldmap_zpad;
for dim = 1:2
    img_fieldmap_interp = 1 / sqrt(size(img_fieldmap_interp,dim)) * fftshift(fft(ifftshift(img_fieldmap_interp, dim), [], dim), dim);
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%%
if read_sign_fieldmap
    img_fieldmap = flip(img_fieldmap,1);
    img_fieldmap_interp = flip(img_fieldmap_interp,1);

    x_fieldmap = flip(x_fieldmap,1);
    y_fieldmap = flip(y_fieldmap,1);
    z_fieldmap = flip(z_fieldmap,1);

    x_fieldmap_interp = flip(x_fieldmap_interp,1);
    y_fieldmap_interp = flip(y_fieldmap_interp,1);
    z_fieldmap_interp = flip(z_fieldmap_interp,1);
end

if phase_sign_fieldmap
    img_fieldmap = flip(img_fieldmap,2);
    img_fieldmap_interp = flip(img_fieldmap_interp,2);

    x_fieldmap = flip(x_fieldmap,2);
    y_fieldmap = flip(y_fieldmap,2);
    z_fieldmap = flip(z_fieldmap,2);

    x_fieldmap_interp = flip(x_fieldmap_interp,2);
    y_fieldmap_interp = flip(y_fieldmap_interp,2);
    z_fieldmap_interp = flip(z_fieldmap_interp,2);
end

%% Display images

slice_number = 24;

figure;
surf(x_fieldmap(:,:,slice_number) * 1e3, y_fieldmap(:,:,slice_number) * 1e3, z_fieldmap(:,:,slice_number) * 1e3, abs(img_fieldmap(:,:,slice_number,1,1,1)), 'EdgeColor', 'none');
axis image;
colormap(gray(256));
xlim([-260 260]);
ylim([-260 260] / 2);
%zlim([-260 260]);
view(0, 90);

figure;
surf(x_fieldmap_interp(:,:,slice_number) * 1e3, y_fieldmap_interp(:,:,slice_number) * 1e3, z_fieldmap_interp(:,:,slice_number) * 1e3, abs(img_fieldmap_interp(:,:,slice_number,1,1,1)), 'EdgeColor', 'none');
axis image;
colormap(gray(256));
xlim([-260 260]);
ylim([-260 260] / 2);
view(0, 90);

%%
close all;
figure;
imagesc(abs(img_fieldmap(:,:,slice_number,1,1,1)));
axis image;

figure;
imagesc(abs(img_fieldmap_interp(:,:,slice_number,1,1,1)));
axis image;




return

%% Perform phase unwrapping
tstart = tic; fprintf('%s: Performing phase unwrapping... ', datetime);
img_phase = angle(img_fieldmap_interp); % [rad]
img_phase = unwrap(img_phase, [], 6); % 2 * Nkx x Nky x nr_slices x 1 x 1 x 2
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Normalize the phase
tstart = tic; fprintf('%s: Normalize the unwrapped phase of multi-echo images... ', datetime);
img_phase = img_phase - repmat(img_phase(:,:,:,:,:,1), [1 1 1 1 1 2]);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform a linear least-squares fit of the phase-time curve
%--------------------------------------------------------------------------
% y in [rad], a in [rad/sec], TE in [sec]
% y1 = a * TE1 + b    [y1]   [TE1 1] [a]
% y2 = a * TE2 + b => [y2] = [TE2 1] [b] => y = A * x
% yn = a * TEn + b    [yn]   [TEn 1]
%
% x_ls = inv(A.' * A) * A.' * y
%--------------------------------------------------------------------------
yy = reshape(permute(img_phase, [6 1 2 3 4 5]), [2 Nkx * Nky * nr_slices]); % [rad]
A = cat(2, TE, ones(2,1));
tstart = tic; fprintf('%s: Performing linear least-squares fitting... ', datetime);
x_ls = inv(A.' * A) * (A.' * yy); % 2 x Nkx * Nky * nr_slices
fieldmap_raw = reshape(x_ls(1,:) / (2 * pi), [Nkx Nky nr_slices]); % [Hz]
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate a threshold mask
img_abs = sum(abs(img_fieldmap_interp(:,:,:,1,1,1)),6);
mask = (img_abs > mean(img_abs(:)) * 0.60);
mask = bwareaopen(mask, 60); % Keep only blobs with an area of 60 pixels or more.

%% Apply a threshold mask
fieldmap_interp_masked = fieldmap_raw .* mask;

%% Calculate the regularization term for the B0 field inhomogeneity
weights = zeros(Nkx, Nky, 'single');
for idx2 = 1:Nky
    for idx1 = 1:Nkx
        %------------------------------------------------------------------
        % Calculate the k-space weight for B0 field inhomogeneity
        %------------------------------------------------------------------
        kx = (-floor(Nkx/2) + idx1 - 1) / Nkx;
        ky = (-floor(Nky/2) + idx2 - 1) / Nky;
        weights(idx1,idx2) = 1 / (1 + w_fB0 * (kx^2 + ky^2))^h_fB0;
    end
end

%% Calculate k-space data
%--------------------------------------------------------------------------
% Siemens: k-space <=> image space
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Applying inverse FFT to move from image space to k-space... ', datetime);
ksp = fieldmap_interp_masked;
for dim = 1:2
    ksp = sqrt(size(ksp,dim)) * fftshift(ifft(ifftshift(ksp, dim), [], dim), dim);
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform smoothing via k-space filtering
ksp = bsxfun(@times, weights, ksp);

%% Calculate a smoothed fieldmap [Hz]
%--------------------------------------------------------------------------
% Siemens: k-space <=> image space
%--------------------------------------------------------------------------
tstart = tic; fprintf('%s: Applying forward FFT to move from k-space to image space... ', datetime);
fieldmap_smooth = ksp;
for dim = 1:2
    fieldmap_smooth = 1 / sqrt(size(fieldmap_smooth,dim)) * fftshift(fft(ifftshift(fieldmap_smooth, dim), [], dim), dim);
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
fieldmap_smooth = real(fieldmap_smooth);

%% Display a raw fieldmap
figure('Color', 'w');
montage(fieldmap_raw, 'DisplayRange', []);
colormap(colorcet('D1'));
colorbar;
clim([-15 15]);
export_fig(fullfile(output_path, 'fieldmap_montage_raw'), '-r300', '-tif');

%% Display a masked, interpolated fieldmap
figure('Color', 'w');
montage(fieldmap_interp_masked, 'DisplayRange', []);
colormap(colorcet('D1'));
colorbar;
clim([-15 15]);
export_fig(fullfile(output_path, 'fieldmap_montage_masked'), '-r300', '-tif');

%% Display a filtered, masked, interpolated fieldmap
figure('Color', 'w');
montage(fieldmap_smooth, 'DisplayRange', []);
colormap(colorcet('D1'));
colorbar;
clim([-15 15]);
export_fig(fullfile(output_path, 'fieldmap_montage_masked_filtered'), '-r300', '-tif');

%% Write a .cfl file
for slice_number = 1:nr_slices
    %----------------------------------------------------------------------
    % fieldmap_smooth (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('fieldmap_smooth_slc%d_gnc%d_%s', slice_number, gnc_flag, slice_type);
    cfl_file = fullfile(output_path, img_filename);
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, fieldmap_smooth(:,:,slice_number));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end
