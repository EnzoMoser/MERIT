[points, axes_] = merit.domain.get_pix_xys(20, 0.08);

    frequencies = dlmread('example_data/frequencies.csv');
    antenna_locations = dlmread('example_data/antenna_locations.csv');
    channel_names = dlmread('example_data/channel_names.csv');
    
    scan1 = dlmread('example_data/B0_P3_p000.csv'); %#ok<*DLMRD>
    scan2 = dlmread('example_data/B0_P3_p036.csv');
    
    signals = scan1-scan2;


delays = merit.beamform.get_delays(channel_names, antenna_locations(:, 1:2), ...
  relative_permittivity=8);
% Perform imaging
img = abs(merit.beamform(signals, frequencies, points, delays, ...
        merit.beamformers.DAS));

% Convert to grid for image display
grid_ = merit.domain.img2grid(img, points, axes_{:});

%im_slice = merit.visualize.get_slice(img, points, axes_);
figure;
imagesc(axes_{:}, grid_)


