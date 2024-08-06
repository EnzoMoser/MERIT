function pixel_delay_time = get_delays(prop_speed, m_size, pix_xs, pix_ys, antenna_locations)
% NxMXJ
% Returns an array of the time it takes for each signal to reach each 2D
% point.

%% Calculate the delay for each signal at each position
number_antennas = size(antenna_locations, 1);

% Get the ant x/y positions
ant_xs = antenna_locations(:, 2);
ant_ys = antenna_locations(:, 1);

% Init array for storing pixel time-delays
p_ds = zeros([number_antennas, m_size, m_size]);

for a_pos = 1:number_antennas  % For each antenna position

    % Find x/y position differences of each pixel from antenna
    x_diffs = pix_xs - ant_xs(a_pos);
    y_diffs = pix_ys - ant_ys(a_pos);
    
    % Calculate one-way time-delay of propagation from antenna to
    % each pixel
    p_ds(a_pos, :, :) = sqrt(x_diffs.^2 + y_diffs.^2);
end

pixel_delay_time = p_ds ./ prop_speed; % Convert distance delay to time delay.

% Apply extra time delay for monostatic. Constant taken from Reimer.
pixel_delay_time = pixel_delay_time + 0.19e-9;
pixel_delay_time = -pixel_delay_time;
end
