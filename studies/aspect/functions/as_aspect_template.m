function [template, psi, b2] = as_aspect_template(az_deg, height_m, uv, aspect)
%AS_ASPECT_TEMPLATE Bird aspect factor of each gate (dev).
%   [TEMPLATE, PSI, B2] = AS_ASPECT_TEMPLATE(AZ_DEG, HEIGHT_M, UV, ASPECT)
%   returns the factor by which the radar overestimates the density of a
%   gate seen at azimuth AZ_DEG and height HEIGHT_M (m above the antenna):
%       log T = b2 (-cos 2 psi) + b1 cos psi - log N,  psi = az - track(h),
%   with the track (direction of motion) of the volume's velocity profile
%   UV ([n_layer x 2] u, v of 100 m layers above the antenna) at the
%   gate's height, so that birds seen side-on (psi = +-90 deg) are
%   brighter (b2 > 0) and tail-on (psi = 0, flying away from the radar)
%   dimmer than head-on (b1 < 0). N makes the azimuthal mean of T one: a
%   uniform field keeps its azimuthal mean density, so the bird
%   cross-section used for the conversion keeps its meaning (an average
%   over viewing directions).
%   The side-on amplitude grows with height and speed (overlap fit):
%       b2 = max(0, ASPECT.b2 + ASPECT.b2_per_km (h - ASPECT.h0_m) / 1000
%                   + ASPECT.b2_per_ms (|uv(h)| - ASPECT.v0_ms)),
%   b1 = ASPECT.b1. Gates in layers without a velocity get T = 1.
%   AZ_DEG and HEIGHT_M have the same size; TEMPLATE, PSI (deg) and B2 too.

n_layer = size(uv, 1);
layer = min(n_layer, max(1, floor(height_m / 100) + 1));
% Layers without a profile value take the nearest layer with one
valid = all(isfinite(uv), 2) & hypot(uv(:, 1), uv(:, 2)) > 0;
template = ones(size(az_deg));
psi = nan(size(az_deg));
b2 = zeros(size(az_deg));
if ~any(valid)
    return
end
valid_id = find(valid);
nearest = interp1(valid_id, valid_id, (1:n_layer)', "nearest", "extrap");
uv_filled = uv(nearest, :);
track = atan2d(uv_filled(:, 1), uv_filled(:, 2));
speed = hypot(uv_filled(:, 1), uv_filled(:, 2));
layer_height = ((1:n_layer)' - 0.5) * 100;
b2_layer = max(0, aspect.b2 + aspect.b2_per_km * (layer_height - aspect.h0_m) / 1000 + ...
    aspect.b2_per_ms * (speed - aspect.v0_ms));
% Front/back term, optionally ramped from zero at the antenna to full at b1_ramp_m
b1_layer = aspect.b1 * ones(n_layer, 1);
if isfield(aspect, "b1_ramp_m") && aspect.b1_ramp_m > 0
    b1_layer = aspect.b1 * min(1, layer_height / aspect.b1_ramp_m);
end

% Normalisation of each layer: azimuthal mean of exp(log T) = 1
psi_grid = (0.5:1:359.5);
log_norm = log(mean(exp(b2_layer .* -cosd(2 * psi_grid) + b1_layer .* cosd(psi_grid)), 2));

psi = mod(az_deg - track(layer), 360);
b2 = b2_layer(layer);
template = exp(b2 .* -cosd(2 * psi) + b1_layer(layer) .* cosd(psi) - log_norm(layer));
end
