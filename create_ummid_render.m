clear
clc

m_size = 100;

%% ADI Radius - For some reason not in the metadata
ADI_RADS.A1 = 3.87e-2;
ADI_RADS.A2 = 5e-2;
ADI_RADS.A3 = 5.66e-2;
ADI_RADS.A11 = 4.48e-2;
ADI_RADS.A12 = 4.57e-2;
ADI_RADS.A13 = 4.84e-2;
ADI_RADS.A14 = 5.19e-2;
ADI_RADS.A15 = 5.53e-2;
ADI_RADS.A16 = 5.74e-2;

%% Range of Information
roi_rad = 8e-2; % The radius of information is 8cm.

%% Choose scan number to work on
scan_num = 1;
%% Load data
% Metadata
metadata = load('../um_bmid/datasets/gen-three/clean/md_list_s11_adi.mat');
% Extract the metadata array from the struct dytpe it is stored in
metadataFieldName = fieldnames(metadata);
metadata = getfield(metadata, metadataFieldName{1}); %#ok<GFLD>

% Signal data
scan_data = load('../um_bmid/datasets/gen-three/clean/fd_data_s11_adi.mat');
% Extract the data array from the struct dtype it is stored in
scanDataFieldName = fieldnames(scan_data);
scan_data = getfield(scan_data, scanDataFieldName{1}); %#ok<GFLD>

%% Setup frequencies
frequencies = linspace( 1e9, 9e9, size(scan_data, 2) ); % Frequencies
frequencies = frequencies(:);

%% Shrink signal frequencies - for better performance
% Shrink the frequencies for scan_data
frequency_ids = frequencies >=2e9; % Only include frequencies above this number.
frequencies = frequencies(frequency_ids);
% Down-sample freq. Only use every [sample_divide]-th element.
sample_divide = 12; % Number to divide the sample by.
frequencies = frequencies(1:sample_divide:end);
% Make sure signal frequencies match
scan_data = scan_data(:, frequency_ids, :);
scan_data = scan_data(:, 1:sample_divide:end, :);

%% Load One Signal
% Get the full signal
org_signal = squeeze(scan_data(scan_num, :, :));
org_metadata = metadata{scan_num};

%% Load metadata from signal
% Get antenna phase radius
ant_rad = metadata{scan_num}.ant_rad * 1e-2; % The distance from the antenna to the centre.
ant_phase_rad = ant_rad + 2.4e-2; % Add the distance from the centre to the antenna's phase centre. "2.4cm" was measured by Tyson Reimer. 
% Get breast tumor
if isnan(org_metadata.tum_diam)
    disp('This scan does NOT have a tumor!');
    tum_coords = [0, 0];
    tum_rad = 0;
else
    tum_coords = [ org_metadata.tum_x, org_metadata.tum_y ] * 1e-2;
    tum_rad = ( org_metadata.tum_diam / 2 ) * 1e-2;
end
% Get breast adipose
adi_coords = [ org_metadata.adi_x, org_metadata.adi_y ] * 1e-2;
split_adi_id = strsplit(org_metadata.phant_id, 'F'); % Split the string after every "F". This will tell us the adi ID
adi_rad = ADI_RADS.(split_adi_id{1}); % Radius of the adi-pose. For some reason, it is not included in the metadata.

% %% Plot the signal
% figure;
% data_channel1 = [ org_signal(:, 1) ] ;
% channel1_magnitude = mag2db(abs(data_channel1));
% channel1_phase = unwrap(angle(data_channel1));
% subplot(2, 1, 1);
% plot(frequencies, channel1_magnitude);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% legend('clean signal');
% subplot(2, 1, 2);
% plot(frequencies, channel1_phase);
% xlabel('Frequency (Hz)');
% ylabel('Phase (rad)');
% legend('clean signal');

%% Calculate signal speed
prop_speed = um_bmid.get_signal_speed(...
    ant_phase_rad, ...
    adi_rad, ...
    adi_coords=adi_coords, ...
    tum_coords=tum_coords, ...
    tum_rad=tum_rad, ...
    m_size=m_size ...
    );

%% Calculate the antenna locations
number_antennas = size(scan_data, 3); % The number of antenna locations.
starting_antenna_angle = deg2rad(-130);

antenna_locations = um_bmid.get_antenna_locations(...
    ant_phase_rad, ...
    number_antennas=number_antennas, ...
    starting_antenna_angle=starting_antenna_angle...
    );

%% Create meshgrid for calculations
%Define the x/y points on each axis
xs = linspace(-roi_rad, roi_rad, m_size);
ys = linspace(roi_rad, -roi_rad, m_size);
% Cast these to 2D for ease
[ pix_xs, pix_ys ] = meshgrid(xs, ys);

% Find the distance from each pixel to the center of the model space
pix_dist_from_center = sqrt(pix_xs.^2 + pix_ys.^2);

%% Calculate region of interest
% Get the region of interest as all the pixels inside the
% circle-of-interest
in_roi = false([m_size, m_size]);
in_roi(pix_dist_from_center < roi_rad) = true;

%% Calculate time delay for each signal at each point
pixel_delay_time = um_bmid.get_delays(prop_speed, m_size, pix_xs, pix_ys, antenna_locations);

%% Calculate phase factor
pix_ts = pixel_delay_time;
phase_fac = exp(-2i * pi * frequencies(:) .* reshape(pix_ts, 1, size(pix_ts, 1), size(pix_ts, 2), size(pix_ts, 3)));

%% Perform DAS beamform
% Convert adi_cal_cropped to have singleton dimensions for broadcasting
adi_cal_cropped_expanded = org_signal;
% Perform the element-wise multiplication and summation using implicit expansion
das_adi_recon_temp1 = adi_cal_cropped_expanded .* phase_fac.^2;
das_adi_recon = sum(sum(das_adi_recon_temp1, 1), 2);
% Squeeze the result to remove singleton dimensions
das_adi_recon = squeeze(das_adi_recon);

% Apply extra time delay for monostatic. Constant taken from Reimer.
delays = merit.beamform.get_delays([1:72; 1:72]', antenna_locations, relative_permittivity=( 299792458 ./ prop_speed ).^2 );

apple = @(ps) delays(ps) - 2*(0.19e-9);

% Generate imaging domain
[points, axes_] = merit.domain.get_pix_xys(m_size, roi_rad);

new_pix_ts = apple(points);

new_das = merit.beamform(org_signal, frequencies(:), points, apple, um_bmid.DAS);

% Convert to grid for image display
grid_ = merit.domain.img2grid(squeeze(abs(new_das)), points, axes_{:});

%% Plot image
img = rot90(grid_, 2);
img_to_plt = img .* ones(size(img), "like", img);
% Set the pixels outside the antenna trajectory to NaN
img_to_plt( ~in_roi ) = nan;

figure()  % Make the figure window
hold on

axes_{1} = xs;
axes_{2} = ys;

imagesc(axes_{:}, img_to_plt);
colormap("jet");
colorbar
title("merit")
axis equal
lims = clim;
hold off

%% Plot image
img = squeeze(abs(das_adi_recon));
img_to_plt = img .* ones(size(img), "like", img);
% Set the pixels outside the antenna trajectory to NaN
img_to_plt( ~in_roi ) = nan;

figure()  % Make the figure window
hold on

axes_{1} = xs;
axes_{2} = ys;

imagesc(axes_{:}, img_to_plt);
colormap("jet");
colorbar
title("reimer")
axis equal
clim(lims)
hold off