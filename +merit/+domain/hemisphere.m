function [points, axes_] = hemisphere(options)
  % [points, axes_] = merit.domain.hemisphere('radius', r, 'resolution', res)
  %   points is a list of points in a hemisphere with resolution spacing between them.
  %   axes_ is the discrete points in each direction.
  %   r is the radius of the hemisphere
  %   resolution is the spacing between points and can be different for each axis.
  arguments
      options.resolution (1,1){mustBeNumeric, mustBeReal} = 1e-3                        %resolution must be numeric, real value, default 1e-3
      options.radius (1,1){mustBeNumeric, mustBeReal, mustBeScalarOrEmpty} = 7e-2       %radius must be a numeric, real scalal, default value 7e-2
      options.x {numeric_or_rvector}  = []                                    % x, y, z must be numeric or real vectors
      options.y {numeric_or_rvector}  = []                    
      options.z {numeric_or_rvector}  = []
      options.no_z = false                                                              % no_z (no z axis) is default false
  end


  resolution = options.resolution;
  if isscalar(resolution)
    resolution = repmat(resolution, [1, 3]);
  end

  radius_ = options.radius+5e-3;

  full_axis = @(r) -radius_:r:radius_;
  half_axis = @(r) 0:r:radius_;

  if options.no_z == false
      axes_ = {options.x, options.y, options.z};
      if isempty(axes_{1})
        axes_{1} = full_axis(resolution(1));
      end
      if isempty(axes_{2})
        axes_{2} = full_axis(resolution(2));
      end
      if isempty(axes_{3})
        axes_{3} = half_axis(resolution(3));
      end

      [Xs, Ys, Zs] = ndgrid(axes_{:});
      area_ = Xs.^2 + Ys.^2+Zs.^2 <= options.radius.^2;
  else
      axes_ = {options.x, options.y};
      if isempty(axes_{1})
        axes_{1} = full_axis(resolution(1));
      end
      if isempty(axes_{2})
        axes_{2} = full_axis(resolution(2));
      end

      [Xs, Ys] = ndgrid(axes_{:});
      area_ = Xs.^2 + Ys.^2 <= options.radius.^2;
  end

  points = merit.beamform.imaging_domain(area_, axes_{:});
end

function numeric_or_rvector(a)
    if ~isempty(a)
        validateattributes(a, {'numeric'}, {'vector', 'real'})
    end
end

