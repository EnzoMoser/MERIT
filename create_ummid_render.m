clear
clc

m_size = 150;

%% Choose scan number to work on
scan_num = 1;

ADI_RADS.A1 = 3.87e-2;
ADI_RADS.A2 = 5e-2;
ADI_RADS.A3 = 5.66e-2;
ADI_RADS.A11 = 4.48e-2;
ADI_RADS.A12 = 4.57e-2;
ADI_RADS.A13 = 4.84e-2;
ADI_RADS.A14 = 5.19e-2;
ADI_RADS.A15 = 5.53e-2;
ADI_RADS.A16 = 5.74e-2;

%% Load data
metadata = load('../um_bmid/datasets/gen-three/clean/md_list_s11_adi.mat');
% Extract the metadata array from the struct dytpe it is stored in
metadataFieldName = fieldnames(metadata);
metadata = getfield(metadata, metadataFieldName{1}); %#ok<GFLD>
org_metadata = metadata{scan_num};

scan_data = load('../um_bmid/datasets/gen-three/clean/fd_data_s11_adi.mat');
% Extract the data array from the struct dtype it is stored in
scanDataFieldName = fieldnames(scan_data);
scan_data = getfield(scan_data, scanDataFieldName{1}); %#ok<GFLD>

%% Cut extra data from scan_data
scan_data = scan_data(1:8, :, :); % Limit it to 7 scans and 2 antennas

%% Setup frequencies
frequencies = linspace( 1e9, 9e9, size(scan_data, 2) ); % Frequencies
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

if isnan(org_metadata.tum_diam)
    disp('This scan does NOT have a tumor!');
end

% Get the adipose-only version.
    cropped_signal =org_signal;

% Get the adipose-plus-fiboglandular-only version (everything except tumor)
%adifib_signal = squeeze(scan_data(org_metadata.fib_ref_id, :, :));

%% Plot the acquired scans.
figure;
data_channel1 = [ cropped_signal(:, 1) ] ;
channel1_magnitude = mag2db(abs(data_channel1));
channel1_phase = unwrap(angle(data_channel1));
subplot(2, 1, 1);
plot(frequencies, channel1_magnitude);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('clean signal');
%title(sprintf('Channel (%d, %d) Magnitude', channel_names(1, :)));
subplot(2, 1, 2);
plot(frequencies, channel1_phase);
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
legend('clean signal');
%title(sprintf('Channel (%d, %d) Phase', channel_names(1, :)));
% 

%% Get signal radius
ant_rad = metadata{scan_num}.ant_rad * 1e-2; % The distance from the antenna to the centre.
sig_rad = ant_rad + 2.4e-2; % "2.4cm" is by Reimer. The distance from the centre to the antenna's phase centre.

split_adi_id = strsplit(org_metadata.phant_id, 'F'); % Split the string after every "F".
adi_rad = ADI_RADS.(split_adi_id{1}); % Radius of the adi-pose. For some reason, it is not included in the metadata.

roi_rad = 8e-2; % The radius of information is 8cm.

%Define the x/y points on each axis
full_xs = linspace(-sig_rad, sig_rad, m_size);
full_ys = linspace(sig_rad, -sig_rad, m_size);

% Cast these to 2D for ease later
[ full_pix_xs, full_pix_ys ] = meshgrid(full_xs, full_ys);

% Find the distance from each pixel to the center of the model space
pix_dist_from_center = sqrt(full_pix_xs.^2 + full_pix_ys.^2);

% Get the region of interest as all the pixels inside the
% circle-of-interest
full_in_roi = false([m_size, m_size]);
full_in_roi(pix_dist_from_center < sig_rad) = true;

% % This is needed to calculate the permitivity and hence speed, but for
% % now we're skipiing this part.
air_permitivity = 1;
adi_permitivity = 7.08;
tum_permitivity = 77.11;
skin_permitivity = 40;
skin_thickness = 0;

tum_x = org_metadata.tum_x * 1e-2;
tum_rad = org_metadata.tum_diam * 1e-2 / 2;
tum_y = org_metadata.tum_y * 1e-2;
adi_x = org_metadata.adi_x * 1e-2;
adi_y = org_metadata.adi_y * 1e-2;

% Compute the pixel distances from the center of each tissue
% component (excluding skin)
pix_dist_from_adi = sqrt((full_pix_xs - adi_x).^2 + (full_pix_ys - adi_y).^2);
pix_dist_from_tum = sqrt((full_pix_xs - tum_x).^2 + (full_pix_ys - tum_y).^2);

% Initialize breast model to be uniform background medium
breast_model = air_permitivity .* ones([m_size, m_size]);

% Assign the tissue components
breast_model(pix_dist_from_adi < (adi_rad + skin_thickness) ) = skin_permitivity;
breast_model(pix_dist_from_adi < adi_rad) = adi_permitivity;
breast_model(pix_dist_from_tum < tum_rad) = tum_permitivity;

vac_speed = 3e8;

prop_speed = mean( ( vac_speed ./ sqrt( breast_model(full_in_roi) ) ), "all" );

number_antennas = size(scan_data, 3); % The number of antenna locations.

starting_antenna_angle = deg2rad(-130);
antenna_angles = (linspace(0, (1 - (1/number_antennas) ) * 2 * pi, number_antennas)); % If number_antennas = 72, then steps of 5 from 0 to 355 (but in radians)
antenna_angles = antenna_angles + starting_antenna_angle; % The starting antenna angle is not always 0 degrees.

antenna_locations = permute ( [ ( cos(antenna_angles) * sig_rad ); ( sin(antenna_angles) * sig_rad ) ], [2,1] );

%Define the x/y points on each axis
xs = linspace(-roi_rad, roi_rad, m_size);
ys = linspace(roi_rad, -roi_rad, m_size);

% Cast these to 2D for ease latersig_rad
[ pix_xs, pix_ys ] = meshgrid(xs, ys);
% Find the distance from each pixel to the center of the model space
pix_dist_from_center = sqrt(pix_xs.^2 + pix_ys.^2);

% Get the region of interest as all the pixels inside the
% circle-of-interest
in_roi = false([m_size, m_size]);
in_roi(pix_dist_from_center < roi_rad) = true;

% Get the ant x/y positions
ant_xs = antenna_locations(:, 2);
ant_ys = antenna_locations(:, 1);

% Init array for storing pixel time-delays
p_ds = zeros([number_antennas, m_size, m_size]);

for a_pos = 1:number_antennas  % For each antenna position

    % Find x/y position differences of each pixel from antenna
    x_diffs = pix_xs - ant_xs(a_pos);
    y_diffs = pix_ys - ant_ys(a_pos);
    
    % Calculate one-way time-delay of propagation from antenna to
    % each pixel
    p_ds(a_pos, :, :) = sqrt(x_diffs.^2 + y_diffs.^2);
end

pixel_delay_time = p_ds ./ prop_speed; % Convert distance delay to time delay, relative to permittivity.

% Apply extra time delay for monostatic. Constant taken from Reimer.
pixel_delay_time = pixel_delay_time + 0.19e-9;

% Assuming recon_fs is a vector of size [N, 1]
% and pix_ts is a 3D array of size [72, 150, 150]

recon_fs = frequencies(:);
pix_ts = pixel_delay_time;

phase_fac = exp(-1i * 2 * pi * recon_fs .* reshape(pix_ts, 1, size(pix_ts, 1), size(pix_ts, 2), size(pix_ts, 3)));

% Convert adi_cal_cropped to have singleton dimensions for broadcasting
adi_cal_cropped_expanded = reshape(cropped_signal, size(cropped_signal, 1), size(cropped_signal, 2), 1, 1);

% Perform the element-wise multiplication and summation using implicit expansion
das_adi_recon = sum(sum(adi_cal_cropped_expanded .* phase_fac.^(-2), 1), 2);

% Squeeze the result to remove singleton dimensions
das_adi_recon = squeeze(das_adi_recon);

img = abs(das_adi_recon);

img_to_plt = img .* ones(size(img), "like", img);

% Set the pixels outside the antenna trajectory to NaN
img_to_plt( ~in_roi ) = nan;

% Bounds for x/y axes ticks in plt
tick_bounds = [-roi_rad, roi_rad, -roi_rad, roi_rad];

figure()  % Make the figure window
hold on

axes_{1} = xs;
axes_{2} = ys;

imagesc(axes_{:}, img_to_plt);
colormap("jet");
colorbar
axis equal
hold off


% channel_one = 1:number_antennas;
% %channel_two = mod ( channel_one+11, number_antennas) + 1; % The receiving antenna is 60 degrees ahead of the transmitting antenna
% channel_names = permute ( [ channel_one; channel_one ], [2, 1] );
% 
% % Calculate delays
% % merit.get_delays returns a function that calculates the delay
% %   to each point from every antenna.
% 
% delays = merit.beamform.get_delays(channel_names, antenna_locations, ...
%   'relative_permittivity', av_perm);
% 
% % Perform imaging
% beamformer = merit.beamformers.DAS;
% img = abs(merit.beamform(org_signal, frequencies, points, delays, beamformer));

% Convert to grid for image display
%grid_ = merit.domain.img2grid(img, points, axes_{:});

%figure
%imagesc(axes_{:}, grid_);
%axis equal