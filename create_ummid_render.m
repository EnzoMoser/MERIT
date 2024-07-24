clear
clc

%% Choose scan number to work on
scan_num = 1;

ADI_RADS.A1 = 3.87;
ADI_RADS.A2 = 5;
ADI_RADS.A3 = 5.66;
ADI_RADS.A11 = 4.48;
ADI_RADS.A12 = 4.57;
ADI_RADS.A13 = 4.84;
ADI_RADS.A14 = 5.19;
ADI_RADS.A15 = 5.53;
ADI_RADS.A16 = 5.74;

%% Load data
metadata = load('../um_bmid/datasets/gen-three/simple-clean/matlab-data/metadata_gen_three.mat');
% Extract the metadata array from the struct dytpe it is stored in
metadataFieldName = fieldnames(metadata);
metadata = getfield(metadata, metadataFieldName{1}); %#ok<GFLD>
org_metadata = metadata{scan_num};

if isnan(org_metadata.tum_diam)
    error('This scan does not have a tumor! Select a different scan!');
end

scan_data = load('../um_bmid/datasets/gen-three/simple-clean/matlab-data/fd_data_gen_three_s11.mat');
% Extract the data array from the struct dtype it is stored in
scanDataFieldName = fieldnames(scan_data);
scan_data = getfield(scan_data, scanDataFieldName{1}); %#ok<GFLD>

%% Setup frequencies
frequencies = linspace( 1e9, 8e9, size(scan_data, 2) ); % Frequencies
% Shrink the frequencies for scan_data
frequency_ids = frequencies >=2e9; % Only include frequencies above this number.
frequencies = frequencies(frequency_ids);
% Down-sample freq. Only use every [sample_divide]-th element.
sample_divide = 12; % Number to divide the sample by.
frequencies = frequencies(:, 1:12:end, :);

% Make sure signal frequencies match
scan_data = scan_data(:, frequency_ids, :);
scan_data = scan_data(:, 1:sample_divide:end, :);

%% Setup signals
% Get the full signal
org_signal = squeeze(scan_data(scan_num, :, :));

% Get the adipose-only version.
adi_signal = squeeze(scan_data(org_metadata.adi_ref_id, :, :));

% Get the adipose-plus-fiboglandular-only version (everything except tumor)
%adifib_signal = squeeze(scan_data(org_metadata.fib_ref_id, :, :));

%% Plot the acquired scans.
figure;
data_channel1 = [ org_signal(:, 1), ( org_signal(:, 1) - adi_signal(:, 1) )];
channel1_magnitude = mag2db(abs(data_channel1));
channel1_phase = unwrap(angle(data_channel1));
subplot(2, 1, 1);
plot(frequencies, channel1_magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('signal', 'clean signal');
%title(sprintf('Channel (%d, %d) Magnitude', channel_names(1, :)));
subplot(2, 1, 2);
plot(frequencies, channel1_phase);
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
legend('signal', 'clean signal');
%title(sprintf('Channel (%d, %d) Phase', channel_names(1, :)));
% 

%% Get signal radius
ant_rad = metadata{scan_num}.ant_rad * 1e-2; % The distance from the antenna to the centre.
sig_rad = ant_rad + 2.4e-2; % "2.4cm" is by Reimer. The distance from the centre to the antenna's phase centre.

split_adi_id = strsplit(org_metadata.phant_id, 'F'); % Split the string after every "F".
adi_rad = ADI_RADS.(split_adi_id{1}); % Radius of the adi-pose. For some reason, it is not included in the metadata.

% % This is needed to calculate the permitivity and hence speed, but for
% % now we're skipiing this part.
% air_permitivity = 1;
% adi_permitivity = 7.08;
% fib_permitivity = 44.94;
% tum_permitivity = 77.11;
% skin_permitivity = 40;
% skin_thickness = 0;
% 
% pixel_dist_from_adi = [];
% pixel_size = 0;
% breast_model = air_permitivity * ones(pixel_size);
% breast_model[ (pixel_dist_from_adi < adi_rad) ] = skin_perm;

prop_speed = 288739680.6368426; % Ripped it from the Python calculation

number_antennas = size(scan_data, 3); % The number of antenna locations.

starting_antenna_angle = deg2rad(-130);
antenna_angles = (linspace(0, (1 - (1/number_antennas) ) * 2 * pi, number_antennas)); % If number_antennas = 72, then steps of 5 from 0 to 355 (but in radians)
antenna_angles = antenna_angles + starting_antenna_angle; % The starting antenna angle is not always 0 degrees.

antenna_locations = permute ( [ ( cos(antenna_angles) * ant_rad ); ( sin(antenna_angles) * ant_rad ) ], [2,1] );

roi_rad = 8e-2; % The radius of information is 8cm.
% Generate imaging domain
[points, axes_] = merit.domain.hemisphere('radius', roi_rad, 'resolution', 2e-3, 'no_z', true);
pixel_size = size(points, 1);
pixel_delay_dis = zeros(number_antennas, pixel_size);

for a_loc = 1:number_antennas

    x_diff = (points(:, 1) - squeeze( antenna_locations(a_loc, 1) ) );
    y_diff = (points(:, 2) - squeeze( antenna_locations(a_loc, 2) )  );

    pix_dis = sqrt( x_diff.^2 + y_diff.^2 );

    pixel_delay_dis(a_loc, :) = pix_dis;

end

pixel_delay_time = pixel_delay_dis ./ prop_speed;

% Apply extra time delay for monostatic. Constant taken from Reimer.
pixel_delay_time = pixel_delay_time + 0.19e-9;

phase_factor = exp(-2i * pi * pixel_delay_time);

cropped_signal = org_signal - adi_signal;

p_size = size(points, 1);

das_reconstruction = zeros([p_size 1]);
for point = 1:p_size
    delayed_signal = zeros(size(cropped_signal));
    for i = 1:size(phase_factor, 1)
        delayed_signal(:, i) = cropped_signal(:, i) .* phase_factor(i, point);
    end
    das_reconstruction(point) = shiftdim(sum(sum(delayed_signal, 2).^2, 1), 2);
end

%channel_one = 1:number_antennas;
%channel_two = mod ( channel_one+11, number_antennas) + 1; % The receiving antenna is 60 degrees ahead of the transmitting antenna
%channel_names = permute ( [ channel_one; channel_two ], [2, 1] );



% % Calculate delays
% % merit.get_delays returns a function that calculates the delay
% %   to each point from every antenna.
% 
% delays = merit.beamform.get_delays(channel_names, antenna_locations, ...
%   'relative_permittivity', 8);
% 
% % Perform imaging
% beamformer = merit.beamformers.DAS;
% img = abs(merit.beamform(org_signal, frequencies, points, delays, beamformer));

% Convert to grid for image display
grid_ = merit.domain.img2grid(abs(das_reconstruction), points, axes_{:});

figure
imagesc(axes_{:}, grid_);
axis equal