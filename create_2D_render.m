% Not accurate, as it's squashing 3D data into 2D.

if ~exist('frequencies', 'var') % Check if frequencies exist.
    frequencies = dlmread('example_data/frequencies.csv');
    antenna_locations = dlmread('example_data/antenna_locations.csv');
    channel_names = dlmread('example_data/channel_names.csv');
    
    scan1 = dlmread('example_data/B0_P3_p000.csv'); %#ok<*DLMRD>
    scan2 = dlmread('example_data/B0_P3_p036.csv');
    
    signals = scan1-scan2;
end

% Generate imaging domain
[points, axes_] = merit.domain.hemisphere('radius', 7e-2, 'resolution', 2.5e-3, 'no_z', true);

% Calculate delays
% merit.get_delays returns a function that calculates the delay
%   to each point from every antenna.
delays = merit.beamform.get_delays(channel_names, antenna_locations(:, 1:2), ...
  'relative_permittivity', 8);

% Perform imaging
img = abs(merit.beamform(signals, frequencies, points, delays, ...
        merit.beamformers.DAS));

% Convert to grid for image display
grid_ = merit.domain.img2grid(img, points, axes_{:});

%im_slice = merit.visualize.get_slice(img, points, axes_);
figure
imagesc(axes_{:}, grid_);