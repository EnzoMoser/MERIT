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

number_of_scans = size(metadata, 2);

relative_permittivity = zeros(number_of_scans, 1);

for scan_num = 1:number_of_scans
    % Load one signal
    org_metadata = metadata{scan_num};

    % Load metadata from signal
    %--- Get antenna phase radius
    ant_rad = metadata{scan_num}.ant_rad * 1e-2; % The distance from the antenna to the centre.
    ant_phase_rad = ant_rad + 2.4e-2; % Add the distance from the centre to the antenna's phase centre. "2.4cm" was measured by Tyson Reimer. 
    
    %--- Get the breast tumor coords and radius
    if isnan(org_metadata.tum_diam)
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
    
    % Calculate signal speed
    prop_speed = um_bmid.get_signal_speed(...
        ant_phase_rad, ...
        adi_rad, ...
        adi_coords=adi_coords, ...
        tum_coords=tum_coords, ...
        tum_rad=tum_rad ...
        );
    
    % Calculate delays
    c_0 = 299792458; % Vaccuum speed, taken from "merit.beamform.get_delays.m"
    relative_permittivity(scan_num)=( c_0 ./ prop_speed ).^2;
    disp("Relative Permittivity MEAN:")
    disp(mean(relative_permittivity))
end