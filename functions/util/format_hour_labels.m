function labels = format_hour_labels(time_values)

labels = compose("t=%dh", round(hours(time_values)));
end
