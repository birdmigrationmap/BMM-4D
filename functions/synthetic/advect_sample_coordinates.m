function [x_adv, y_adv] = advect_sample_coordinates(x, y, H, gz, uv, dt_days, reference_id)
%ADVECT_SAMPLE_COORDINATES Translate samples into a common advection frame.
% reference_id = 0 moves scans to a frame one step before the first scan.
% A positive scan index moves every scan into that scan's time frame and is
% the preferred choice when covariance and simulation should share one
% physical reference time.

if nargin < 7
    reference_id = ceil(size(uv, 2) / 2);
end

x_adv = x;
y_adv = y;
dt_seconds = dt_days(:)' * 24 * 60 * 60;

for i_t = 1:size(uv, 2)
    if reference_id == 0
        step_ids = 1:i_t;
        sign_step = -1;
    elseif i_t < reference_id
        step_ids = i_t:reference_id-1;
        sign_step = 1;
    else
        step_ids = reference_id:i_t-1;
        sign_step = -1;
    end

    dx = zeros(size(H(:, i_t)));
    dy = zeros(size(H(:, i_t)));

    for i_step = step_ids
        uv_sample = interp1(gz(:), uv(:, i_step), H(:, i_t), "linear");
        dx = dx + sign_step .* real(uv_sample) .* dt_seconds(i_step);
        dy = dy + sign_step .* imag(uv_sample) .* dt_seconds(i_step);
    end

    x_adv(:, i_t) = x(:, i_t) + dx;
    y_adv(:, i_t) = y(:, i_t) + dy;
end
