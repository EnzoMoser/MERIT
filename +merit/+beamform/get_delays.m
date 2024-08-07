function [calculate_time] = get_delays(channels, antennas, relative_permittivity)
arguments
  % list of channels
  channels {mustBeInteger}
  % list of antenna locations
  antennas {mustBeNumeric}
  % relative permittivity must be a nummeric scaler >= 1
  relative_permittivity {mustBeGreaterThanOrEqual(relative_permittivity,1)}
end
  c_0 = 299792458;

  speed = c_0./sqrt(relative_permittivity);

  antennas = antennas';
  
  function [time] = calculate_(pointsf)
    points = permute(pointsf, [2, 3, 1]);
    distances = sqrt(sum( (antennas - points).^2, 1) );
    time = - ( distances(:, channels(:, 1), :) + distances(:, channels(:, 2), :) ) / speed;
  end
  calculate_time = @calculate_;
end
