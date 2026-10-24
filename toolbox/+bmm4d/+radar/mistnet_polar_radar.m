function [preds, classes, probs] = mistnet_polar_radar(radar)
%MISTNET_POLAR Classify sample volumes in polar coordinates using mistnet

% MatConvNet built with -ffinite-math-only (wsrlib maca64 builds before
% the fix) returns 0 for max pooling of negative values, which silently
% changes MistNet's predictions.
persistent pooling_checked
if isempty(pooling_checked)
    if exist("vl_nnpool", "file") && vl_nnpool(single([-1 -2; -3 -4]), [2 2]) ~= -1
        error("mistnet_polar_radar:BrokenPooling", ...
            "MatConvNet's vl_nnpool returns wrong maxima for negative values. Rebuild wsrlib's MatConvNet (vl_compilenn) with -ffast-math -fno-finite-math-only.")
    end
    pooling_checked = true;
end

% Run the upstream MistNet implementation on its Cartesian g.
[PREDS, PROBS, classes, y, x, elevs] = mistnet(radar);

% Map Cartesian predictions back to the original polar sample volumes.
F = griddedInterpolant({y, x, elevs}, PREDS, 'nearest', 'none');

preds = cell(numel(radar.dz.sweeps),1);
probs = cell(numel(radar.dz.sweeps),1);
for i_elev = 1:numel(radar.dz.sweeps)
    % Convert sample volume coordinate matrices to XYZ
    [~, range, az] = sweep2mat(radar.dz.sweeps(i_elev));
    [RANGE, AZ, ELEV] = ndgrid(range, az, radar.dz.sweeps(i_elev).elev);
    [X, Y, ~] = radar2xyz(RANGE, AZ, ELEV);
    ELEV = min(max(ELEV, elevs(1)), elevs(end));

    % Interpolate predictions onto polar volume
    preds{i_elev} = F(Y, X, ELEV);

    % If class probabilities are requested, also interpolate those, one at a time
    if nargout > 2
        n_classes = size(PROBS, 4);
        probs{i_elev} = nan([size(ELEV) n_classes]);
        for c = 1:n_classes
            F = griddedInterpolant({y, x, elevs}, PROBS(:,:,:,c), 'nearest', 'none');
            probs{i_elev}(:,:,c) = F(Y, X, ELEV);
        end
    end
end
end
