function uv_profiles_interval = average_velocity_time_profile(uv_profiles_time)

uv_profiles_interval = 0.5 .* (uv_profiles_time(:, :, 1:end-1) + uv_profiles_time(:, :, 2:end));
end
