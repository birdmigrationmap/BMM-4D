function fig = create_diagnostic_figure(width_px, height_px, varargin)
% Create a larger fixed-size figure for diagnostics and exploratory exports.

scale_factor = 1.5;
fig = create_export_figure(round(scale_factor * width_px), round(scale_factor * height_px), varargin{:});
end
