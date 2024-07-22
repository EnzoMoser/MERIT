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