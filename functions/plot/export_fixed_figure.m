function export_fixed_figure(fig, output_file, varargin)
% Export rasters at the figure canvas size; delegate vector/append exports.

resolution = 300;
content_type = "";
append_export = false;
background_color = [];

for i_arg = 1:2:numel(varargin)
    name = string(varargin{i_arg});
    value = varargin{i_arg + 1};
    if name == "Resolution"
        resolution = value;
    elseif name == "ContentType"
        content_type = string(value);
    elseif name == "Append"
        append_export = value;
    elseif name == "BackgroundColor"
        background_color = value;
    end
end

[~, ~, ext] = fileparts(output_file);
if append_export || lower(string(ext)) == ".pdf" || content_type == "vector"
    exportgraphics(fig, output_file, varargin{:})
    return
end

old_color = fig.Color;
old_invert = fig.InvertHardcopy;
old_paper_units = fig.PaperUnits;
old_paper_position = fig.PaperPosition;
old_paper_size = fig.PaperSize;
cleanup = onCleanup(@() set(fig, ...
    "Color", old_color, ...
    "InvertHardcopy", old_invert, ...
    "PaperUnits", old_paper_units, ...
    "PaperPosition", old_paper_position, ...
    "PaperSize", old_paper_size));
if ~isempty(background_color)
    fig.Color = background_color;
end
fig.InvertHardcopy = "off";
fig.PaperUnits = "inches";
fig.PaperPosition = [0 0 fig.Position(3:4)];
fig.PaperSize = fig.Position(3:4);

if any(lower(string(ext)) == [".jpg" ".jpeg"])
    device = "-djpeg";
else
    device = "-dpng";
end

print(fig, output_file, device, "-r" + resolution)
normalize_raster_size(output_file, round(fig.Position(3) * resolution), round(fig.Position(4) * resolution), fig.Color)
end

function normalize_raster_size(output_file, target_width, target_height, background_color)

image_data = imread(output_file);
if size(image_data, 2) == target_width && size(image_data, 1) == target_height
    return
end

common_height = min(target_height, size(image_data, 1));
common_width = min(target_width, size(image_data, 2));

if ismatrix(image_data)
    fill_value = cast(round(255 * mean(background_color)), class(image_data));
    canvas = repmat(fill_value, target_height, target_width);
    canvas(1:common_height, 1:common_width) = image_data(1:common_height, 1:common_width);
else
    canvas = zeros(target_height, target_width, size(image_data, 3), class(image_data));
    fill_color = cast(round(255 * background_color), class(image_data));
    for i_channel = 1:size(image_data, 3)
        canvas(:, :, i_channel) = fill_color(i_channel);
    end
    canvas(1:common_height, 1:common_width, :) = image_data(1:common_height, 1:common_width, :);
end

imwrite(canvas, output_file)
end
