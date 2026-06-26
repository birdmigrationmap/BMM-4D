function uv_profiles_time = sample_synthetic_velocity_profile(grid, time_offsets, varargin)

p = inputParser;
addRequired(p, "grid");
addRequired(p, "time_offsets");
addParameter(p, "speed_min_kmh", 30);
addParameter(p, "speed_max_kmh", 50);
parse(p, grid, time_offsets, varargin{:});

n_time = numel(time_offsets);
nz = grid.sz(3);
speed_kmh = p.Results.speed_min_kmh + (p.Results.speed_max_kmh - p.Results.speed_min_kmh) .* rand;
direction_rad = 2 .* pi .* rand;
uv_profile_ms = (speed_kmh ./ 3.6) .* [cos(direction_rad) sin(direction_rad)];

uv_profiles_time = repmat(reshape(uv_profile_ms, [1 2 1]), nz, 1, n_time);
end
