function y = vl_nnpool(x, pool, varargin)
%VL_NNPOOL MATLAB fallback for MatConvNet pooling on Apple Silicon.

if nargin >= 3 && isnumeric(varargin{1})
    error("vl_nnpool:BackwardNotSupported", ...
        "This MATLAB fallback only supports forward pooling for inference.")
end

if isscalar(pool)
    pool = [pool pool];
end

pad = [0 0 0 0];
stride = [1 1];
method = "max";

i = 1;
while i <= numel(varargin)
    key = lower(string(varargin{i}));
    switch key
        case "pad"
            pad = varargin{i + 1};
            i = i + 2;
        case "stride"
            stride = varargin{i + 1};
            i = i + 2;
        case "method"
            method = lower(string(varargin{i + 1}));
            i = i + 2;
        case {"cudnn", "nocudnn"}
            i = i + 1;
        otherwise
            error("vl_nnpool:UnsupportedOption", "Unsupported option: %s", varargin{i})
    end
end

if isscalar(stride)
    stride = [stride stride];
end
if isscalar(pad)
    pad = repmat(pad, 1, 4);
end

orig_size = size(x);
if ndims(x) < 4
    x = reshape(x, [orig_size ones(1, 4 - ndims(x))]);
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
    otherwise
        error("vl_nnpool:UnsupportedMethod", "Unsupported pooling method: %s", method)
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

switch numel(orig_size)
    case 2
        y = reshape(y, size(y, 1), size(y, 2));
    case 3
        y = reshape(y, size(y, 1), size(y, 2), size(y, 3));
end
end
