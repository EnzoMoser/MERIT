
if ~exist('grid_', 'var') % Only create the grid if it doesn't already exist.
frequencies = dlmread('example_data/frequencies.csv');
antenna_locations = dlmread('example_data/antenna_locations.csv');
channel_names = dlmread('example_data/channel_names.csv');

scan1 = dlmread('example_data/B0_P3_p000.csv'); %#ok<*DLMRD>
scan2 = dlmread('example_data/B0_P3_p036.csv');

signals = scan1-scan2;

[points, axes_] = merit.domain.hemisphere(resolution=2.5e-3, radius=7e-2);

%% Calculate delays
% merit.get_delays returns a function that calculates the delay
%   to each point from every antenna.
delays = merit.beamform.get_delays(channel_names, antenna_locations, ...
  relative_permittivity=8);

%% Perform imaging

img = abs(merit.beamform(signals, frequencies, points, delays, ...
        merit.beamformers.DAS));

%% Convert to grid for image display
grid_ = merit.domain.img2grid(img, points, axes_{:});

end

% Check if figure 4 exists
if ishandle(4)
    %clf(4); % Clear figure 4
end
merit.visualize.display_3D_scan(...
    grid_...
    )
