function [colormap_data, c_lim] = colormap_wct(name, n_points)

% Set default n_points if not provided
if nargin == 0
    disp("Possible values: 'dz', 'vr', 'rh', 'ph', 'dr', 'sw'.")
    return
end

if nargin < 2
    n_points = 256;
end

switch name
    case "dz"
        colors = [
            75   235  235  235  255;
            70   153  85   201  255;
            65   255  0    255  255;
            60   192  0    0    255;
            55   214  0    0    255;
            50   255  0    0    255;
            45   255  144  0    255;
            40   231  192  0    255;
            35   255  255  0    255;
            30   0    144  0    255;
            25   0    200  0    255;
            20   0    255  0    255;
            15   0    0    246  255;
            10   1    160  246  255;
            5    0    236  236  255;
            0    187  255  255  255;
            -5   174  238  238  255;
            -10  150  205  205  255;
            -15  102  139  139  255;
            -20  50   79   79   255;
            -25   0   0    0    255
            ];
    case "vr"
        colors = [
            70   255 0   0;
            60   255 0   0;
            50   210 0   0;
            40   180 0   0;
            30   150 0   0;
            20   115 0   0;
            10   70  40  40;
            0    70  70  70;
            -10  40  70  40;
            -20  0   115 0;
            -30  0   150 0;
            -40  0   180 0;
            -50  0   210 0;
            -60  0   255 0;
            -70  0   255 0
            ];
    case "sw"
        colors = [
            30   255 255    0;
            27   231 183    0;
            24   208 112    0;
            21   231  56    0;
            18   255   0    0;
            15   127  93    0;
            12    0 187    0;
            9    78  171   78;
            6    156  156  156;
            3    137  137  137;
            0    118  118  118
            ];

    case "dr"
        colors = [
            20.0  0 100 0;
            7.9  50 205 50;
            7.9  255 255 255;
            7.0  255 255 255;
            6.0  255 255 255;
            5.0  240 120 180;
            4.0  175 0 0;
            3.0  220 10 5;
            2.0  255 255 100;
            1.5  90 221 98;
            1.0  68 248 212;
            0.25 10 10 155;
            0.0  142 121 181;
            0.0  220 220 220;
            -2.0 165 165 165;
            -4.0 110 110 110;
            -6.0 55 55 55;
            -7.9 0 0 0;
            -7.9 185 156 107;
            -13.0 133 87 35
            ];
    case "ph"
        colors = [
            180.0  255 140 255;
            120.0  90 0 0;
            90.0   255 125 0;
            75.0   255 255 0;
            60.0   30 100 0;
            45.0   0 255 0;
            30.0   10 20 95;
            15.0   210 210 180;
            0.0    255 255 255
            ];
    case "rh"
        colors = [
            105 255 180 215;
            100 255 180 215;
            100 139 30 77;
            97 225 3 0;
            95 255 140 0;
            90 255 255 0;
            85 135 215 10;
            80 95 245 100;
            75 120 120 255;
            60 10 10 190;
            45 15 15 140;
            20 100 100 100;
            ];
end

[colormap_data, c_lim] = convert_colormap_wct(colors, n_points);

end

function [colormap_data, c_lim] = convert_colormap_wct(colors, n_points)

% Normalize RGB values
rgb_values = colors(:, 2:4) / 255;
values = colors(:, 1) + 100*eps.*(1:numel(colors(:, 1)))';  % First column for interpolation

c_lim = [min(values)  max(values)];

% Generate colormap indices based on the first column (values)
xq = linspace(c_lim(1), c_lim(2), n_points)';  % Query points for interpolation

% Interpolate for each RGB channel based on the 'values' column
colormap_data = [
    interp1(values, rgb_values(:, 1), xq, 'linear'), ... % Red
    interp1(values, rgb_values(:, 2), xq, 'linear'), ... % Green
    interp1(values, rgb_values(:, 3), xq, 'linear')    % Blue
    ];
end


