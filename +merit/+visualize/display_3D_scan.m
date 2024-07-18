function display_3D_scan(grid_, options)
arguments
    grid_ (:, :, :)

    % Specify figure to display
    options.figure_number (1, 1) {mustBeNumeric} ...
        = 0 % Default value is don't specify.

    % Shift the opacity of the upper quartile between 0 and 1.
    % The pixels will/should never be fully opaque, as that just looks bad.
    options.upper_opacity (1, 1) {mustBeNumeric, mustBePositive, ...
        mustBeLessThanOrEqual(options.upper_opacity, 1)} ...
        = 0.4 % Default value

    % Change the factor for the opacity of the lower quartile by a...
    % factor of the upper_opactiy
    % lower_opacity = upper_opacity * lower_opacity_factor
    options.lower_opacity_factor (1, 1) {mustBeNumeric, mustBePositive, ...
        mustBeLessThanOrEqual(options.lower_opacity_factor, 1)} ...
        = 0.1 % Default value
end

if options.figure_number < 1
    figure;
else
    figure(options.figure_number);
end
hold on;

[lx, ly, lz] = size(grid_);

% Generate the grid
[x, y, z] = meshgrid(1:lx, 1:ly, 1:lz);

% Plot the ball using slice

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

vector_grid = grid_(grid_ ~= 0 & ~isnan(grid_));

max_value = max(vector_grid);
min_value = min(vector_grid);

Q1 = quantile(vector_grid, 0.25);
Q3 = quantile(vector_grid, 0.75);

% Determine the values at specific steps
st = min_value;
md1 = Q1;
md2 = Q3;
ed = max_value;

sep2 = options.upper_opacity;

sep1 = sep2 * options.lower_opacity_factor;
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