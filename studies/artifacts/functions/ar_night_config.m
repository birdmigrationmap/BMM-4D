function c = ar_night_config(night_name)
%AR_NIGHT_CONFIG Settings of a full us_east night of the night set (dev).
%   C = AR_NIGHT_CONFIG(NIGHT_NAME) returns the settings ns_process_nights
%   cleaned the night with (us_east radars, grid, filter, solve grid and
%   domain; production cleaning without the shape correction), for the
%   row NIGHT_NAME of studies/nights/figures/NS_A_nights.csv.

project_root = setup_paths();
nights = readtable(fullfile(project_root, "studies", "nights", "figures", "NS_A_nights.csv"), TextType="string");
row = nights(nights.name == night_name, :);
assert(height(row) == 1, "ar_night_config:Night", "%s is not in NS_A_nights.csv.", night_name)
times = datetime(row.first_frame):minutes(row.step_min):datetime(row.last_frame);
us = case_settings("us_east");
c = bmm4d.config(night_name, "radars", us.radars, "times", times, "grid_bounds", us.grid_bounds, ...
    "filter", us.filter, "solve_aggregate", us.solve_aggregate, "domain", us.domain, ...
    "read_options", {"shape", false});
end
