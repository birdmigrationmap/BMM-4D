%% Sandbox: column steps where the ground crosses a layer boundary
% The field is layered in height above sea level (100 m). A cell is active
% when its top is above the ground (build_active_mask), and the lowest cell
% is about as dense as the one above, so a column that counts whole cells
% jumps by a layer of birds where the ground crosses a layer boundary:
% blocky column maps that follow the terrain (README, "Review of the two
% nights"). Counting the lowest cell by its share above the ground
% (bmm4d.grid.air_fraction) makes the column continuous.
% Also tested, and not kept (no measurable change): the same share in the
% forward operator (birds below the ground do not echo), and the prior's
% radar profiles fitted and applied in height above the ground (raw DEM or
% smoothed over 30 km); results in the README.
% Score: mean log2 change of the column between neighbouring columns
% (4 km) where the ground is one layer higher, and without a step.
% Cases: the sandbox (21 Oct 2026, 6 frames) and the 12 May 2021 night.
% Outputs: figures/ground_<case>.png (sandbox_case_figure, columns with the
% share), ground.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
case_names = ["sandbox" "night_20210512"];

%% Solve and score --------------------------------------------------------
scores = table();
for case_name = case_names
    S_case = sandbox_load(case_name);
    result = sandbox_run(S_case, struct("name", "production " + case_name, "video", false));
    cell_air = bmm4d.grid.air_fraction(S_case.g, S_case.active_id);
    whole = ground_step(S_case, result.field, ones(size(cell_air)));
    share = ground_step(S_case, result.field, cell_air);
    scores = [scores; table(case_name, whole(1), whole(2), share(1), share(2), VariableNames=["case" ...
        "whole_cells_no_step" "whole_cells_one_layer" "air_share_no_step" "air_share_one_layer"])]; %#ok<AGROW>
    sandbox_case_figure(S_case, {result}, fullfile(figure_dir, "ground_" + case_name + ".png"), ...
        ceil(size(result.field, 2) / 2), "Columns with the air share of the lowest cell: " + case_name);
end
disp(scores)
writetable(scores, fullfile(figure_dir, "ground.csv"))

%% Local functions ----------------------------------------------------------
function step = ground_step(S, field, cell_weight)
% Mean log2 column change between x-neighbours (4 km), mean over frames:
% [no ground step, ground one layer higher]; columns above 20 birds km^-2,
% each cell counted by CELL_WEIGHT.
g = S.g;
[ix, iy, iz] = ind2sub(g.sz, S.active_id);
column_id = sub2ind(g.sz(1:2), ix, iy);
lowest = reshape(accumarray(column_id, iz, [prod(g.sz(1:2)) 1], @min, NaN), g.sz(1:2));
d_lowest = lowest(2:end, :) - lowest(1:end - 1, :);
n_time = size(field, 2);
step = zeros(1, 2);
for i_t = 1:n_time
    column = reshape(accumarray(column_id, double(field(:, i_t)) .* cell_weight, [prod(g.sz(1:2)) 1]) * g.dz / 1000, g.sz(1:2));
    d_column = log2(max(column(2:end, :), 1) ./ max(column(1:end - 1, :), 1));
    ok = ~isnan(d_lowest) & column(2:end, :) > 20 & column(1:end - 1, :) > 20;
    step = step + [mean(d_column(ok & d_lowest == 0)), mean(d_column(ok & d_lowest == 1))] / n_time;
end
end
