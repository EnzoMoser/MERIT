function [calculate_time] = get_delays(channels, antennas, options)
arguments
  channels                                                                                                                      % lsit of channels
  antennas                                                                                                                          % list of antenna locations
  options.relative_permittivity {mustBeNumeric, mustBeScalarOrEmpty, mustBeGreaterThanOrEqual(options.relative_permittivity,1)}     % relative permittivity must be a nummeric scaler >= 1
end
  c_0 = 299792458;
  relative_permittivity = options.relative_permittivity;

  speed = c_0./sqrt(relative_permittivity);

  antennas = antennas';
  
  function [time] = calculate_(pointsf)
    points = permute(pointsf, [2, 3, 1]);
    distances = sqrt(sum( (antennas - points).^2, 1) );
    time = - ( distances(:, channels(:, 1), :) + distances(:, channels(:, 2), :) ) / speed;
  end
  calculate_time = @calculate_;
end
