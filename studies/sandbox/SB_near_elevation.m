%% Sandbox: do the elevations agree near the radar?
% In the reconstructions the density within 15 km of a radar is 25-30%
% below that at 15-60 km at 500-1500 m: dark cones above the radars, with
% bright spikes along single beams, and on 12 May dense towers in the
% lowest kilometre (README, "Review of the two nights"). Test on the data:
% solve without the rows within 30 km of their radar (row weight 0), then
% compare those held-out rows with what the solution predicts for them,
% by range and elevation. The prediction comes from the farther beams of
% the same radar and from the neighbours, at the same heights.
% If the ratio depends on range and not on elevation, it is the gate size
% (the gate estimator's range bias, studies/gate_density GD_I); if it
% depends on elevation at the same range, it is the beam or the clutter.
% Second view: the same ratio for rows 3-30 km against the expected number
% of birds per gate, lambda = predicted density x counting volume
% (bmm4d.radar.gate.beam), by elevation group: if the gate estimator is
% biased when gates hold few birds, the groups fall on one curve.
% Cases: the sandbox (21 Oct 2026, 6 frames) and the 12 May 2021 night.
% Outputs: figures/near_elevation.png, near_elevation.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
case_names = ["sandbox" "night_20210512"];
range_edges_km = [2 4 6 8 10 12.5 15 20 25 30 40 60];
elevation_edges = [0 0.7 1.1 1.5 2.1 2.8 3.5 4.5 5.5 7];
elevation_label = ["0.5°" "0.9°" "1.3°" "1.8°" "2.4°" "3.1°" "4.0°" "5.1°" "6.4°"];

%% Held-out near rows -----------------------------------------------------------
summary = table();
ratios = cell(numel(case_names), 1);
beam = bmm4d.radar.gate.beam();
lambda_edges = 10 .^ (-1:0.25:2.5);
lambda_centre = sqrt(lambda_edges(1:end - 1) .* lambda_edges(2:end));
group_edges = [0 1.5 2.8 7];
group_label = ["< 1.5°" "1.5-2.8°" "≥ 2.8°"];
lambda_ratios = cell(numel(case_names), 1);
for i_c = 1:numel(case_names)
    S_case = sandbox_load(case_names(i_c));
    result = sandbox_run(S_case, struct("name", "rows from 30 km " + case_names(i_c), ...
        "row_weight", @(r, e, z) double(r >= 30000), "video", false));
    obs = S_case.observations;
    predicted = obs.F(:, S_case.active_id) * double(result.field);
    observed = full(obs.values);
    has = observed ~= 0 & predicted > 1;
    range_bin = discretize(double(obs.row_range_m) / 1000, range_edges_km);
    elevation_bin = discretize(double(obs.row_elev_deg), elevation_edges);
    ratio = nan(numel(range_edges_km) - 1, numel(elevation_label));
    count = zeros(size(ratio));
    for i_r = 1:numel(range_edges_km) - 1
        for i_e = 1:numel(elevation_label)
            rows = range_bin == i_r & elevation_bin == i_e;
            use = has(rows, :);
            o = observed(rows, :);
            p = predicted(rows, :);
            count(i_r, i_e) = nnz(use);
            if count(i_r, i_e) >= 200
                ratio(i_r, i_e) = sum(o(use)) / sum(p(use));
            end
        end
    end
    ratios{i_c} = ratio;
    % By birds per gate, rows 3-30 km
    lambda = predicted .* beam.volume_sample_km3(double(obs.row_range_m));
    lambda_bin = discretize(lambda, lambda_edges);
    group = discretize(double(obs.row_elev_deg), group_edges);
    near = obs.row_range_m >= 3000 & obs.row_range_m < 30000;
    lambda_ratio = nan(numel(lambda_centre), numel(group_label));
    for i_l = 1:numel(lambda_centre)
        for i_g = 1:numel(group_label)
            use = has & near & group == i_g & lambda_bin == i_l;
            if nnz(use) > 300
                lambda_ratio(i_l, i_g) = sum(observed(use)) / sum(predicted(use));
            end
        end
    end
    lambda_ratios{i_c} = lambda_ratio;
    centre = (range_edges_km(1:end - 1) + range_edges_km(2:end))' / 2;
    summary = [summary; [table(repmat(case_names(i_c), numel(centre), 1), centre, VariableNames=["case" "range_km"]) ...
        array2table(round(ratio, 2), VariableNames=cellstr(elevation_label))]]; %#ok<AGROW>
end
disp(summary)
writetable(summary, fullfile(figure_dir, "near_elevation.csv"))

%% Plot results -----------------------------------------------------------
% One panel per case: observed / predicted against range, one line per
% elevation, with the gate study's range bias of the estimator (all scans)
gd = readtable(fullfile(project_root, "studies", "gate_density", "figures", "GD_I_range_bias.csv"));
fig = bmm4d.plot.create_export_figure(3600, 1100, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 2 * numel(case_names), TileSpacing="tight", Padding="tight");
title(t, "Near-range rows held out of the solve: observed / predicted by range and elevation", FontSize=14)
colors = parula(numel(elevation_label) + 1);
for i_c = 1:numel(case_names)
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    centre = (range_edges_km(1:end - 1) + range_edges_km(2:end)) / 2;
    for i_e = 1:numel(elevation_label)
        plot(ax, centre, ratios{i_c}(:, i_e), "-o", Color=colors(i_e, :), LineWidth=1.5, MarkerFaceColor=colors(i_e, :), ...
            DisplayName=elevation_label(i_e))
    end
    plot(ax, gd.range_km, gd.All124Scans_estimator, "k--", LineWidth=2, DisplayName="gate study b(d), estimator")
    yline(ax, 1, "k:", HandleVisibility="off")
    xlim(ax, [0 60]); ylim(ax, [0 2]); set(ax, FontSize=12)
    xlabel(ax, "Range (km)"); ylabel(ax, "Observed / predicted")
    legend(ax, Location="southeast", NumColumns=3, FontSize=9)
    title(ax, case_names(i_c) + ": by range", FontWeight="normal", Interpreter="none")
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    group_colors = colors([1 4 8], :);
    for i_g = 1:numel(group_label)
        plot(ax, lambda_centre, lambda_ratios{i_c}(:, i_g), "-o", Color=group_colors(i_g, :), LineWidth=1.5, ...
            MarkerFaceColor=group_colors(i_g, :), DisplayName=group_label(i_g))
    end
    yline(ax, 1, "k:", HandleVisibility="off")
    set(ax, XScale="log", FontSize=12); xlim(ax, [0.1 300]); ylim(ax, [0 2])
    xlabel(ax, "Expected birds per gate \lambda"); ylabel(ax, "Observed / predicted (3-30 km)")
    legend(ax, Location="southeast", FontSize=10)
    title(ax, case_names(i_c) + ": by birds per gate", FontWeight="normal", Interpreter="none")
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "near_elevation.png"), Resolution=150)
close(fig)
