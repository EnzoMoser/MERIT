function antenna_locations = get_antenna_locations(ant_phase_rad, options)
% Nx2
% Returns an array of coordinates corresponding to each antenna
% locations in 2D space.
arguments
    % Distance from antenna's phase centre to the chamber centre.
    ant_phase_rad (1, 1) {mustBeGreaterThanOrEqual(ant_phase_rad, 0)}
    % These defaults are based off of gen_three defaults from Tyson Reimer.
    % (N). The number of antennas. The antennas are set to be spaced out evenly.
    options.number_antennas (1, 1) {mustBeGreaterThanOrEqual(options.number_antennas, 0)} = 72
    options.starting_antenna_angle = deg2rad(-130)
end

antenna_angles = (linspace(0, (1 - (1/options.number_antennas) ) * 2 * pi, options.number_antennas)); % If number_antennas = 72, then steps of 5 from 0 to 355 (but in radians)
antenna_angles = antenna_angles + options.starting_antenna_angle; % The starting antenna angle is not always 0 degrees.

antenna_locations = permute ( [ ( cos(antenna_angles) * ant_phase_rad ); ( sin(antenna_angles) * ant_phase_rad ) ], [2,1] );
end

