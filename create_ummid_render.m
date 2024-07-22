scan_num = 280;

scan_data = load('../um_bmid/datasets/gen-two/clean/fd_data_s21_adi.mat');
% Extract the data array from the struct dtype it is stored in
scanDataFieldName = fieldnames(scan_data);
scan_data = getfield(scan_data, scanDataFieldName{1}); %#ok<GFLD>

metadata = load('../um_bmid/datasets/gen-two/clean/md_list_s21_adi.mat');
% Extract the metadata array from the struct dytpe it is stored in
metadataFieldName = fieldnames(metadata);
metadata = getfield(metadata, metadataFieldName{1}); %#ok<GFLD>

signals = squeeze(scan_data(scan_num, :, :));
frequencies = linspace( 1e9, 8e9, size(signals, 1) );
radius = metadata{scan_num}.ant_rad * 1e-2; % The radius of the scan.
number_antennas = size(signals, 2); % The number of antenna locations.

antenna_angles = (linspace(0, (355 / 360) * 2 * pi, number_antennas));
antenna_locations = permute ( [ ( cos(antenna_angles) * radius ); ( sin(antenna_angles) * radius ) ], [2,1] );

%% Plot the acquired scans.
figure;
data_channel1 = signals(:, 1);
channel1_magnitude = mag2db(abs(data_channel1));
channel1_phase = unwrap(angle(data_channel1));
subplot(2, 1, 1);
plot(frequencies, channel1_magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Signal');
%title(sprintf('Channel (%d, %d) Magnitude', channel_names(1, :)));
subplot(2, 1, 2);
plot(frequencies, channel1_phase);
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
legend('Signal');
%title(sprintf('Channel (%d, %d) Phase', channel_names(1, :)));

% Generate imaging domain
[points, axes_] = merit.domain.hemisphere('radius', radius, 'resolution', 2.5e-3, 'no_z', true);

% Calculate delays
% merit.get_delays returns a function that calculates the delay
%   to each point from every antenna.

one_channel = 1:number_antennas;
channel_names = permute ( [ one_channel; one_channel ], [2, 1] );

delays = merit.beamform.get_delays(channel_names, antenna_locations, ...
  'relative_permittivity', 8);

% Perform imaging
img = abs(merit.beamform(signals, frequencies, points, delays, ...
        merit.beamformers.DAS));

% Convert to grid for image display
grid_ = merit.domain.img2grid(img, points, axes_{:});

%im_slice = merit.visualize.get_slice(img, points, axes_);
figure
imagesc(axes_{:}, grid_);