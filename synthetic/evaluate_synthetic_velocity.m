function [u, v] = evaluate_synthetic_velocity(X, Y, Z, velocity_model, i_interval)
%EVALUATE_SYNTHETIC_VELOCITY Evaluate a reproducible horizontal velocity field.
% Direction is a compass bearing measured clockwise from north. U is the
% eastward g-x component and V is the northward g-y component.

if velocity_model.model == "uniform"
    sample_grid = X + Y + Z;
    speed = velocity_model.base_speed_mps + zeros(size(sample_grid), "like", sample_grid);
    direction = velocity_model.base_direction_rad;
else
    i_interval = max(1, min(velocity_model.n_interval, i_interval));
    speed = interpolate_velocity_field( ...
        velocity_model.coarse_x, velocity_model.coarse_y, velocity_model.coarse_z, ...
        velocity_model.speed_fields(:, :, :, i_interval), X, Y, Z);
    direction = interpolate_velocity_field( ...
        velocity_model.coarse_x, velocity_model.coarse_y, velocity_model.coarse_z, ...
        velocity_model.direction_fields(:, :, :, i_interval), X, Y, Z);
    speed = max(1, speed);
end

u = speed .* sin(direction);
v = speed .* cos(direction);
end

function field_values = interpolate_velocity_field(coarse_x, coarse_y, coarse_z, field_nodes, X, Y, Z)

[query_x, query_y, query_z] = expand_query_arrays(X, Y, Z);
F = griddedInterpolant({coarse_x(:), coarse_y(:), coarse_z(:)}, field_nodes, "linear", "nearest");
field_values = F(query_x, query_y, query_z);
end

function [query_x, query_y, query_z] = expand_query_arrays(X, Y, Z)

query_template = X + Y + Z;
query_x = zeros(size(query_template), "like", query_template) + X;
query_y = zeros(size(query_template), "like", query_template) + Y;
query_z = zeros(size(query_template), "like", query_template) + Z;
end
