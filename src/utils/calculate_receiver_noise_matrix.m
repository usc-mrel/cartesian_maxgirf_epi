function [Psi,inv_L] = calculate_receiver_noise_matrix(noise_file)


%% Read an ismrmrd file
tic; fprintf('%s: Reading an ISMRMRD file: %s... ', datetime, noise_file);
if exist(noise_file, 'file')
    dset = ismrmrd.Dataset(noise_file, 'dataset');
    fprintf('done! (%6.4f sec)\n', toc);
else
    warning('%s: File %s does not exist.  Please generate it.' , datetime, noise_file);
    Psi = single(1);
    inv_L = single(1);
    return;
end
raw_data = dset.readAcquisition();

%% Sort noise data
acq_is_noise_measurement = raw_data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
meas = raw_data.select(find(acq_is_noise_measurement));
nr_repetitions = length(meas.data); % number of repetitions
[nr_samples,nr_channels] = size(meas.data{1});
noise = complex(zeros(nr_channels, nr_samples, nr_repetitions, 'single'));
for idx = 1:nr_repetitions
    noise(:,:,idx) = meas.data{idx}.'; % nr_samples x nr_channels => nr_channels x nr_samples
end

%% Calculate the receiver noise matrix
% Use the definition in Appendix B of Pruessmann et al. (MRM 46:638–651 (2001))
% ns denotes the number of noise samples taken per channel
% eta lists these samples in an nc x ns matrix
tstart = tic; fprintf('%s: Calculating the receiver noise matrix... ', datetime);
Nc = nr_channels;
Ns = nr_samples * nr_repetitions;
eta = reshape(noise, [Nc Ns]);
Psi = eta * eta' / Ns; % Equation B1
fprintf('done! (%6.4f sec)\n', toc(tstart));

%% Calculate the Cholesky decomposition of the receiver noise matrix
tstart = tic; fprintf('%s: Calculating the Cholesky decomposition of the receiver noise matrix... ', datetime);
L = chol(Psi, 'lower'); % Equation B4
inv_L = inv(L);
fprintf('done! (%6.4f sec)\n', toc(tstart));

end
