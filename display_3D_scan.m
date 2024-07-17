frequencies = dlmread('example_data/frequencies.csv');
antenna_locations = dlmread('example_data/antenna_locations.csv');
channel_names = dlmread('example_data/channel_names.csv');

scan1 = dlmread('example_data/B0_P3_p000.csv');
scan2 = dlmread('example_data/B0_P3_p036.csv');

signals = scan1-scan2;

[points, axes_] = merit.domain.hemisphere('radius', 7e-2, 'resolution', 2.5e-3);

%% Calculate delays
% merit.get_delays returns a function that calculates the delay
%   to each point from every antenna.
delays = merit.beamform.get_delays(channel_names, antenna_locations, ...
  'relative_permittivity', 8);

%% Perform imaging

img = abs(merit.beamform(signals, frequencies, points, delays, ...
        merit.beamformers.DAS));

%% Convert to grid for image display
grid_ = merit.domain.img2grid(img, points, axes_{:});

[lx, ly, lz] = size(grid_);

% Generate the grid
[x, y, z] = meshgrid(1:lx, 1:ly, 1:lz);

% Plot the ball using slice
figure(4);
hold on;

hz = slice(x, y, z, grid_, [], [], 1:lz);
set(hz, 'EdgeColor', 'none');

hy = slice(x, y, z, grid_, [], 1:ly, []);
set(hy, 'EdgeColor', 'none');

hx = slice(x, y, z, grid_, 1:lx, [], []);
set(hx, 'EdgeColor', 'none');

colormap jet; % Set the colormap to 'jet'
colorbar; % Show color bar to indicate density values
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Scan Density Visualization');
view(3); % Set the view to 3D
axis tight; % Fit the axes to the data
grid on; % Turn on the grid
% Add transparency
alpha('color'); % Use color data for transparency

max_value = max(grid_(grid_ ~= 0 & ~isnan(grid_)));
min_value = min(grid_(grid_ ~= 0 & ~isnan(grid_)));

vector_grid = grid_(grid_ ~= 0 & ~isnan(grid_));
Q1 = quantile(vector_grid, 0.25);
Q3 = quantile(vector_grid, 0.75);

% Determine the values at specific steps
st = min_value;
md1 = Q1;
md2 = Q3;
ed = max_value;

sep2 = 0.4;

sep1 = sep2 / 10;
stps1 = ( ed - st) / ( md1 - st );

stps2 = ( ed - st) / ( md2 - md1 );

stps3 = ( ed - st) / ( ed - md2 );

% Adjust alpha_map to reflect the relative differences
alpha_map = [ ...
    linspace( 0.01, sep1, stps1 ) ...
    linspace( sep1, sep2, stps2 )  ...
    linspace( sep2, 1, stps3 ) ...
    ];

% Apply the custom alpha map
alphamap(alpha_map);

axis equal;
hold off;