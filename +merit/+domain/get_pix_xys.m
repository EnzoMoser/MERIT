function[points, axes_] = get_pix_xys(m_size, roi_rho)


%% compute arrays for x values and y values
xs = linspace(-roi_rho, roi_rho, m_size);  % create an array for circle radius roi_rho with resoultion of m_size for x values
ys = -xs;  % create an array for circle radius roi_rho with resoultion of m_size for y values


%% create matrix grid from x and y values
[x_dists, y_dists] = meshgrid(xs, ys); 
axes_ = {xs, ys};

area_ = x_dists.^2 + y_dists.^2 <= roi_rho.^2;

points = merit.beamform.imaging_domain(area_, axes_{:});

end