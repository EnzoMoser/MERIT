scan_num = 45; % 45 and 71 are good tests.

metadata = load('../um_bmid/datasets/gen-three/clean/md_list_s21_adi.mat');
% Extract the metadata array from the struct dytpe it is stored in
metadataFieldName = fieldnames(metadata);
metadata = getfield(metadata, metadataFieldName{1}); %#ok<GFLD>

scan_data = load('../um_bmid/datasets/gen-three/clean/fd_data_s11_adi.mat');
% Extract the data array from the struct dtype it is stored in
scanDataFieldName = fieldnames(scan_data);
scan_data = getfield(scan_data, scanDataFieldName{1}); %#ok<GFLD>

scan1 = squeeze(scan_data(scan_num, :, :));

scan_data = load('../um_bmid/datasets/gen-three/clean/fd_data_s21_adi.mat');
% Extract the data array from the struct dtype it is stored in
scanDataFieldName = fieldnames(scan_data);
scan_data = getfield(scan_data, scanDataFieldName{1}); %#ok<GFLD>

scan2 = squeeze(scan_data(scan_num, :, :));

signals = scan2;

frequencies = linspace( 1e9, 8e9, size(scan1, 1) );
radius = metadata{scan_num}.ant_rad * 1e-2; % The radius of the scan.
number_antennas = size(scan1, 2); % The number of antenna locations.

antenna_angles = (linspace(0, (1 - (1/number_antennas) ) * 2 * pi, number_antennas)); % If number_antennas = 72, then steps of 5 from 0 to 355 (but in radians)
antenna_locations = permute ( [ ( cos(antenna_angles) * radius ); ( sin(antenna_angles) * radius ) ], [2,1] );

channel_one = 1:number_antennas;
channel_two = mod ( channel_one+11, number_antennas) + 1; % The receiving antenna is 60 degrees ahead of the transmitting antenna
channel_names = permute ( [ channel_one; channel_two ], [2, 1] );

%% Plot the acquired scans.
figure;
data_channel1 = [ signals(:, 1)];
channel1_magnitude = mag2db(abs(data_channel1));
channel1_phase = unwrap(angle(data_channel1));
subplot(2, 1, 1);
plot(frequencies, channel1_magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('s21');
title(sprintf('Channel (%d, %d) Magnitude', channel_names(1, :)));
subplot(2, 1, 2);
plot(frequencies, channel1_phase);
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
legend('s21');
title(sprintf('Channel (%d, %d) Phase', channel_names(1, :)));
% 
% Generate imaging domain
[points, axes_] = merit.domain.hemisphere('radius', radius, 'resolution', 2.5e-3, 'no_z', true);

% Calculate delays
% merit.get_delays returns a function that calculates the delay
%   to each point from every antenna.

delays = merit.beamform.get_delays(channel_names, antenna_locations, ...
  'relative_permittivity', 8);

% Perform imaging
beamformer = merit.beamformers.DAS;
img = abs(merit.beamform(signals, frequencies, points, delays, beamformer));

% Convert to grid for image display
grid_ = merit.domain.img2grid(img, points, axes_{:});

figure
imagesc(axes_{:}, grid_);
axis equal