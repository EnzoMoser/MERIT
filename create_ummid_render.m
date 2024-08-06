%% Values that you can quickly change:
scan_num = 1; % Choose scan number to work on (there are 200 in total)
m_size = 100; % Choose the number of pixels for each side. Less is faster:

%% The radius for each adipose - For some reason this is not in the metadata
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

%% Load data
%--- Metadata
% Warning: All values in the metadata are in centimetres!!!
% Make sure to convert them to metres when using.
metadata = load('../um_bmid/datasets/gen-three/clean/md_list_s11_adi.mat');
% Extract the metadata array from the struct dytpe it is stored in
metadataFieldName = fieldnames(metadata);
metadata = getfield(metadata, metadataFieldName{1}); %#ok<GFLD>

%--- Signal data
scan_data = load('../um_bmid/datasets/gen-three/clean/fd_data_s11_adi.mat');
% Extract the data array from the struct dtype it is stored in
scanDataFieldName = fieldnames(scan_data);
scan_data = getfield(scan_data, scanDataFieldName{1}); %#ok<GFLD>

%% Setup frequencies
frequencies = linspace( 1e9, 9e9, size(scan_data, 2) ); % Frequencies
frequencies = frequencies(:);

%% Shrink the number of frequencies - for better performance
%--- Shrink the frequencies for scan_data
frequency_ids = frequencies >=2e9; % Only include frequencies above this number.
frequencies = frequencies(frequency_ids);
%--- Down-sample freq. Only use every [sample_divide]-th element.
sample_divide = 12; % Number to divide the sample by.
frequencies = frequencies(1:sample_divide:end);
%--- Make sure signal frequencies match
scan_data = scan_data(:, frequency_ids, :);
scan_data = scan_data(:, 1:sample_divide:end, :);

%% Load One Signal
%--- Get the full signal
org_signal = squeeze(scan_data(scan_num, :, :));
org_metadata = metadata{scan_num};

%% Load metadata from signal
%--- Get antenna phase radius
ant_rad = metadata{scan_num}.ant_rad * 1e-2; % The distance from the antenna to the centre.
ant_phase_rad = ant_rad + 2.4e-2; % Add the distance from the centre to the antenna's phase centre. "2.4cm" was measured by Tyson Reimer. 

%--- Get the breast tumor coords and radius
if isnan(org_metadata.tum_diam)
    disp('This scan does NOT have a tumor!');
    tum_coords = [0, 0];
    tum_rad = 0;
else
    tum_coords = [ org_metadata.tum_x, org_metadata.tum_y ] * 1e-2;
    tum_rad = ( org_metadata.tum_diam / 2 ) * 1e-2;
end

%--- Get breast adipose
adi_coords = [ org_metadata.adi_x, org_metadata.adi_y ] * 1e-2;
split_adi_id = strsplit(org_metadata.phant_id, 'F'); % Split the string after every "F". This will give us the adipose ID
adi_rad = ADI_RADS.(split_adi_id{1}); % Radius of the adi-pose. For some reason, it is not included in the metadata.

%% Calculate signal speed
prop_speed = um_bmid.get_signal_speed(...
    ant_phase_rad, ...
    adi_rad, ...
    adi_coords=adi_coords, ...
    tum_coords=tum_coords, ...
    tum_rad=tum_rad ...
    );

%% Calculate the antenna locations
number_antennas = size(scan_data, 3); % The number of antenna locations.
starting_antenna_angle = deg2rad(-130);

antenna_locations = um_bmid.get_antenna_locations(...
    ant_phase_rad, ...
    number_antennas=number_antennas, ...
    starting_antenna_angle=starting_antenna_angle...
    );

%% Calculate delays
c_0 = 299792458; % Vaccuum speed, taken from "merit.beamform.get_delays.m"
delays_temp = merit.beamform.get_delays([1:number_antennas; 1:number_antennas]', antenna_locations, relative_permittivity=( c_0 ./ prop_speed ).^2 );

delay_constant = 0.19e-9; % Apply extra time delay for monostatic. Constant taken from Reimer.

% Apply the delay constant. It must be doubled as it is applied for the...
%... distance to send and the distance to receive.
delays = @(points) delays_temp(points) - 2*(delay_constant);

% Generate imaging domain
[points, axes_] = um_bmid.get_pix_xys(m_size, roi_rad);

% Create image
img = abs(merit.beamform(org_signal, frequencies(:), points, delays, um_bmid.DAS));

% For some reason, the image is actually flipped. So we must flip it back to normal:
img = flip(img);

% Convert the image to grid
grid_ = merit.domain.img2grid(img, points, axes_{:});

%% Display image
figure()
imagesc(axes_{:}, grid_);
axis equal
colormap("jet");
colorbar
% Matlab likes to display the Y-Direction backwards. So we tell it not to:
set(gca,'YDir','normal'); set(gca,'XDir','normal');