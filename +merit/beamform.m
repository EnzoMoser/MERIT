function [img, sigs] = beamform(signals, axis_, points, delay, image_, options)
arguments
    signals 
    axis_ 
    points 
    delay 
    image_ 
    options.max_memory {mustBeNumeric} = 1e9
    options.gpu {mustBeNumericOrLogical} = false
    options.window = @(a) a
end
  gpu = options.gpu;
  max_memory = options.max_memory;
  window = options.window;

  if isreal(signals)
    bytes = 8;
  else
    bytes = 16;
  end
  if isa(signals, 'single')
    bytes = bytes/2;
    axis_ = single(axis_);
    points = single(points);
  end

  if ~iscell(delay)
    delay = {delay};
  end

  if gpu
    try
      dev = gpuDevice();
      free_mem = dev.FreeMemory;
      if free_mem > numel(signals)*bytes
        max_memory = dev.FreeMemory/4;
        signals = gpuArray(signals);
        points = gpuArray(points);
      else
        gpu = false;
      end
    catch
      gpu = false;
    end
  end

  max_points = @(mem) floor(mem/numel(signals)/bytes);

  nPoints = size(points, 1);
  points_run = max_points(max_memory);

  img = zeros([nPoints, numel(delay), size_t(signals)], 'like', signals);

  for r = 1:points_run:nPoints
    for d = 1:numel(delay)
      rng = r:min(nPoints, r+points_run-1);
      delayed_signals = merit.process.delay(signals, delay{d}(points(rng, :)), axis_);
      img(rng, d, :) = image_(window(delayed_signals));
    end
  end

  if gpu
    img = gather(img);
    dev = gpuDevice();
    reset(dev);
  end
end

function [s] = size_t(a)
  s = size(a);
  if numel(s) <= 2
    s = [s, 1];
  end
  s = s(3:end);
end
