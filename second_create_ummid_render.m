clc

m_size = 50;

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
scan_data11 = load('../um_bmid/datasets/gen-three/clean/fd_data_s11_adi.mat');
% Extract the data array from the struct dtype it is stored in
scanDataFieldName = fieldnames(scan_data11);
scan_data11 = getfield(scan_data11, scanDataFieldName{1}); %#ok<GFLD>

% Signal data
scan_data21 = load('../um_bmid/datasets/gen-three/clean/fd_data_s21_adi.mat');
% Extract the data array from the struct dtype it is stored in
scanDataFieldName = fieldnames(scan_data21);
scan_data21 = getfield(scan_data21, scanDataFieldName{1}); %#ok<GFLD>

scan_data = scan_data11;

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
% Generate imaging domain
[points, axes_,] = merit.domain.get_pix_xys(m_size, roi_rad);

function [y] = get_delays(prop_speed, m_size, pix_xs, pix_ys, antenna_locations, axes_)

pixel_delay_time = um_bmid.get_delays(prop_speed, m_size, pix_xs, pix_ys, antenna_locations);
xs = axes_{1};
ys = axes_{2};

function [x] = calculate_(list_points)
    positions = zeros(size(list_points));

    for i = 1:size(list_points, 1)
        lx = list_points(i, 1);
        ly = list_points(i, 2);

        % Find the index in all_points that matches (x, y)
        [~, idx] = ismember(lx, xs);
        [~, idy] = ismember(ly, ys);

        if ~isempty(idx) & ~isempty(idy)
            % If found, store the position
            positions(i, :) = [idx, idy];
        else
            error("What")
        end
    end
    positions = permute(positions, [2,1]);

    x = zeros(1, size(pixel_delay_time, 1), size(positions, 2));
    for i = 1:size(positions, 2)
        x(1, :, i) = squeeze(pixel_delay_time(:, positions(1, i), positions(2, i)));
    end
end

y = @calculate_;

end

channels = permute([1:72; 1:72], [2, 1]);

delays = merit.beamform.get_delays(channels, antenna_locations(:, 1:2), ...
  relative_permittivity=1.0802);

% Perform imaging
img = abs(merit.beamform(org_signal, frequencies, points, delays, ...
        merit.beamformers.DAS));

% Convert to grid for image display
grid_ = merit.domain.img2grid(img, points, axes_{:});

figure()
hold on
imagesc(axes_{:}, grid_);

colormap("jet");
colorbar
axis equal
hold off

zz_M = squeeze(delays(points(1:2, :)));

delays = get_delays(prop_speed, m_size, pix_xs, pix_ys, antenna_locations(:, 1:2), axes_);

zz_R = squeeze(delays(points(1:2, :)));

% Perform imaging
img = abs(merit.beamform(org_signal, frequencies, points, delays, ...
        um_bmid.DAS));

% Convert to grid for image display
grid_ = merit.domain.img2grid(img, points, axes_{:});

figure()
hold on
imagesc(axes_{:}, grid_);

colormap("jet");
colorbar
axis equal
hold off