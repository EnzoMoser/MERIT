%antenna_locations = dlmread('example_data/antenna_locations.csv');
%channel_names = dlmread('example_data/channel_names.csv');

scan_data = load('../um_bmid/datasets/gen-two/clean/fd_data_s21_adi.mat');

% Extract the data array from the struct dtype it is stored in
scanDataFieldName = fieldnames(scan_data);
scan_data = getfield(scan_data, scanDataFieldName{1});

signal = squeeze(scan_data(1, :, :));
frequencies = linspace( 1e9, 8e9, size(scan_data, 2) );

%% Plot the acquired scans.
figure;
data_channel1 = signal(:, [1, 10, 20]);
channel1_magnitude = mag2db(abs(data_channel1));
channel1_phase = unwrap(angle(data_channel1));
subplot(2, 1, 1);
plot(frequencies, channel1_magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Original Scan', 'Rotated Scan');
%title(sprintf('Channel (%d, %d) Magnitude', channel_names(1, :)));
subplot(2, 1, 2);
plot(frequencies, channel1_phase);
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
legend('Original Scan', 'Rotated Scan');
%title(sprintf('Channel (%d, %d) Phase', channel_names(1, :)));


% %% Generate imaging domain and visualise
% figure(3)
% [points, axes_] = merit.domain.hemisphere('radius', 7e-2, 'resolution', 2.5e-3);
% subplot(1, 1, 1);
% scatter3(points(:, 1), points(:, 2), points(:, 3), '+');
% 
% %% Calculate delays
% % merit.get_delays returns a function that calculates the delay
% %   to each point from every antenna.
% delays = merit.beamform.get_delays(channel_names, antenna_locations, ...
%   'relative_permittivity', 8);
% 
% %% Perform imaging
% 
% img = abs(merit.beamform(signals, frequencies, points, delays, ...
%         merit.beamformers.DAS));
% 
% %% Convert to grid for image display
% %grid_ = merit.domain.img2grid(img, points, axes_{:});
% 
% im_slice = merit.visualize.get_slice(img, points, axes_, 'z', 35e-3);
% figure(4)
% imagesc(axes_{1:2}, im_slice);