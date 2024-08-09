function merit_data = verify_ummid_render(number_of_scans, m_size, pwd_path, scan_dir, md_dir)
arguments
    number_of_scans {mustBeInteger, mustBePositive}
    m_size {mustBeInteger, mustBePositive}
    pwd_path
    scan_dir
    md_dir
end
% Add the main directory to this scripts path
cd(pwd_path);

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
metadata = load(scan_dir);
% Extract the metadata array from the struct dytpe it is stored in
metadataFieldName = fieldnames(metadata);
metadata = getfield(metadata, metadataFieldName{1}); %#ok<GFLD>

%--- Signal data
scan_data = load(md_dir);
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
%--- Down-sample freq. Only use every "n_divide"-th element.
n_divide = 12; % Number to divide the sample by.
frequencies = frequencies(1:n_divide:end);
%--- Make sure signal frequencies match
scan_data = scan_data(:, frequency_ids, :);
scan_data = scan_data(:, 1:n_divide:end, :);

merit_data = zeros(number_of_scans, m_size, m_size);

for scan_num = 1:number_of_scans
    
    %% Load One Signal
    %--- Get the full signal
    org_signal = squeeze(scan_data(scan_num, :, :));
    org_metadata = metadata{scan_num};
    
    %% Load metadata from signal
    %--- Get antenna phase radius
    % The distance from the antenna to the centre.
    ant_rad = metadata{scan_num}.ant_rad * 1e-2;
    % Add the distance from the centre to the antenna's phase centre. "2.4cm" was measured by Tyson Reimer. 
    ant_phase_rad = ant_rad + 2.4e-2;
    
    %--- Get the breast tumor coords and radius
    if isnan(org_metadata.tum_diam)
        tum_coords = [0, 0];
        tum_rad = 0;
    else
        tum_coords = [ org_metadata.tum_x, org_metadata.tum_y ] * 1e-2;
        tum_rad = ( org_metadata.tum_diam / 2 ) * 1e-2;
    end
    
    %--- Get breast adipose (in metres)
    adi_coords = [ org_metadata.adi_x, org_metadata.adi_y ] * 1e-2;
    % Split the string after every "F". This will give us the adipose ID
    split_adi_id = strsplit(org_metadata.phant_id, 'F');
    % Radius of the adi-pose. For some reason, it is not included in the metadata.
    adi_rad = ADI_RADS.(split_adi_id{1});
    
    %% Calculate the antenna locations
    number_antennas = size(scan_data, 3); % The number of antenna locations.
    starting_antenna_angle = deg2rad(-130);
    
    antenna_locations = merit.domain.create_circumference(...
        ant_phase_rad, ...
        number_antennas, ...
        starting_antenna_angle...
        );
    
    %% Calculate delays
    c_0 = 299792458; % Vaccuum speed, taken from "merit.beamform.get_delays.m"
    % Calculate signal speed
    prop_speed = get_signal_speed(...
        ant_phase_rad, ...
        adi_rad, ...
        adi_coords=adi_coords, ...
        tum_coords=tum_coords, ...
        tum_rad=tum_rad ...
        );
    % Get the permittivity
    relative_permittivity = (c_0 ./ prop_speed).^2;
    
    delays_temp = merit.beamform.get_delays([1:number_antennas; 1:number_antennas]', antenna_locations, relative_permittivity);
    
    %{
    Apply extra time delay for monostatic. Constant taken from T. Reimer's measurements.
    https://github.com/TysonReimer/ORR-EPM/blob/5680df25fae9a3ee0ff3fd0fbb238694efc39a11/umbms/hardware/antenna.py#L41
    %}
    delay_constant = 0.19e-9;
    
    %{
    Apply the delay constant. It must be doubled as it is applied for the
    distance to send and the distance to receive.
    %}
    delays = @(points) delays_temp(points) - 2*(delay_constant);
    
    % Generate imaging domain
    [points, axes_] = merit.domain.get_pix_xys(m_size, roi_rad);
    
    % Create image
    img = abs(merit.beamform(org_signal, frequencies(:), points, delays, DAS));
    
    % For some reason, the image is actually flipped. So we must flip it back to normal:
    img = flip(img);
    
    % Convert the image to grid
    grid_ = merit.domain.img2grid(img, points, axes_{:});
    merit_data(scan_num, :, :) = grid_;
end
end