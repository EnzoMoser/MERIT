function display_3D_scan(grid_, options)
% Display a 3D render of an (:, :, :) grid.
% The density is relative to the highest and lowest numbers.
arguments
    % Matrix containing the data to render.
    grid_ (:, :, :) {mustBeNumeric}

    % --- Optional arguments---

    % Specify figure to display. "0" means don't specify
    options.figure_number (1, 1) {mustBeGreaterThanOrEqual(options.figure_number, 0)} ...
        = 0

    % Change the factor (between 0 and 1) for the max_opactiy of...
    % numbers from the upper quartile to the 90th quartile...
    % ( Default is 1 ).
    options.upper_opacity (1, 1)...
        {mustBeInRange(options.upper_opacity, 0, 1, 'exclude-lower')} ...
        = 1

    % Shift the max_opacity for numbers between the lower and...
    % upper quartile (between 0 and 1).
    % The pixels will/should never be fully opaque, as that just looks bad.
    options.middle_opacity_factor (1, 1) ...% Default value
        {mustBeInRange(options.middle_opacity_factor, 0, 1, 'exclude-lower')} ...
        = 0.4 % Default value

    % Change the factor (between 0 and 1) for the max_opacity of...
    % numbers below the lower quartile by a factor of the upper_opactiy.
    % max_lower_opacity = middle_opacity * lower_opacity_factor
    options.lower_opacity_factor (1, 1) ...
        {mustBeInRange(options.lower_opacity_factor, 0, 1, 'exclude-lower')} ...
        = 0.1
    
    % Specify the colormap. Can take a string or a (:, 3) matrix
    options.colormap_ = "jet"
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

colormap(options.colormap_);
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
Q99 = quantile(vector_grid, 0.90);

% Determine the values at specific steps
st = min_value;
md1 = Q1;
md2 = Q3;
md3 = Q99;
ed = max_value;

sep3 = options.upper_opacity;
sep2 = sep3 * options.middle_opacity_factor;
sep1 = sep2 * options.lower_opacity_factor;
full_diff = ( ed - st);

stps1 = full_diff / ( md1 - st );
stps2 = full_diff / ( md2 - md1 );
stps3 = full_diff / ( md3 - md2 );
stps4 = full_diff / ( ed - md3 );

% Adjust alpha_map to reflect the relative differences
alpha_map = [ ...
    linspace( 0.01, sep1, stps1 ) ...
    linspace( sep1, sep2, stps2 )  ...
    linspace( sep2, sep3, stps3 )  ...
    linspace( sep3, 1, stps4 ) ...
    ];

% Apply the custom alpha map
alphamap(alpha_map);

axis equal;
hold off;