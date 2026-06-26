classdef vlNnPoolFixTest < matlab.unittest.TestCase
    methods (Test)
        function maxPoolingMatchesExplicitLocalMax(testCase)
            addpath(fullfile(fileparts(mfilename("fullpath")), "..", "..", "functions", "wrslib", "matconvnet_fix"), "-begin")
            cleanup = onCleanup(@() rmpath(fullfile(fileparts(mfilename("fullpath")), "..", "..", "functions", "wrslib", "matconvnet_fix"))); %#ok<NASGU>

            rng(1)
            x = single(randn(4, 4, 3));

            y = vl_nnpool(x, [2 2], "pad", [0 1 0 1], "stride", [2 2], "method", "max");
            y_expected = explicitPool(x, [2 2], [0 1 0 1], [2 2], "max");

            testCase.verifyEqual(y, y_expected)
        end

        function avgPoolingMatchesExplicitLocalMean(testCase)
            addpath(fullfile(fileparts(mfilename("fullpath")), "..", "..", "functions", "wrslib", "matconvnet_fix"), "-begin")
            cleanup = onCleanup(@() rmpath(fullfile(fileparts(mfilename("fullpath")), "..", "..", "functions", "wrslib", "matconvnet_fix"))); %#ok<NASGU>

            rng(2)
            x = single(randn(5, 4, 2));

            y = vl_nnpool(x, [3 2], "pad", [1 0 2 1], "stride", [2 1], "method", "avg");
            y_expected = explicitPool(x, [3 2], [1 0 2 1], [2 1], "avg");

            testCase.verifyEqual(y, y_expected, AbsTol=1e-6)
        end
    end
end

function y = explicitPool(x, pool, pad, stride, method)
if ndims(x) < 4
    x = reshape(x, [size(x) ones(1, 4 - ndims(x))]);
end

h = size(x, 1);
w = size(x, 2);
c = size(x, 3);
n = size(x, 4);
out_h = floor((h + pad(1) + pad(2) - pool(1)) / stride(1)) + 1;
out_w = floor((w + pad(3) + pad(4) - pool(2)) / stride(2)) + 1;

switch method
    case "max"
        x_pad = -inf([h + pad(1) + pad(2), w + pad(3) + pad(4), c, n], "like", x);
    case "avg"
        x_pad = zeros([h + pad(1) + pad(2), w + pad(3) + pad(4), c, n], "like", x);
end
x_pad(1 + pad(1):pad(1) + h, 1 + pad(3):pad(3) + w, :, :) = x;

y = zeros(out_h, out_w, c, n, "like", x);
for iy = 1:out_h
    row_idx = (iy - 1) * stride(1) + (1:pool(1));
    for ix = 1:out_w
        col_idx = (ix - 1) * stride(2) + (1:pool(2));
        patch = x_pad(row_idx, col_idx, :, :);
        switch method
            case "max"
                y(iy, ix, :, :) = max(max(patch, [], 1), [], 2);
            case "avg"
                y(iy, ix, :, :) = mean(mean(patch, 1), 2);
        end
    end
end

if size(y, 4) == 1
    y = y(:, :, :, 1);
end
end
