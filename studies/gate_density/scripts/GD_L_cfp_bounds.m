%% Gate density L: the clutter filter's removed power as an upper bound
% On the batch cuts (1.8 deg and above in most scan patterns) the radar's
% clutter filter runs on nearly every gate, and they read about 0.8 of the
% split cuts at the same height (gate-density report, Part III). Putting
% the removed power (CFP) back closes the gap, but also brings back real
% ground clutter. Here a filtered gate is an interval instead: at least the
% output, at most the output plus the removed power (BMM4D.RADAR.GATE.BOUNDS),
% so strong clutter only gives a lower bound. Point-clutter gates (CFP code
% 1, value replaced) are at least the output.
% On the 124 reference scans of GD_J, as GD_I: mean density per scan, sweep
% and range bin for the radar output and three estimators (as in
% production; with the CFP bound; with the CFP bound and point-clutter
% gates), every cut typed split (surveillance pass of a cut scanned twice)
% or batch.
% Results:
% 1. same-height ratio of each cut to the cut below (bin means, the lower
%    cut read where its beam centre is at the same height), by pair type:
%    batch/split should come to split/split and batch/batch if the offset
%    is the filter;
% 2. range bias per group of cuts (GD_FIT_RANGE_BIAS, each group 1 at
%    30-60 km), all scans and scans below 20 birds km^-3, where clutter
%    weighs most against the birds: the lowest cuts near the radar must not
%    rise above 1;
% 3. the share of gate estimates raised more than twofold, by range and cut.
% A fourth estimator puts the bounds on the batch cuts only: the split cuts'
% surveillance pass is filtered only where the radar found clutter, so its
% removed power is mostly clutter (the lowest cuts rose by a third at 30-60
% km with the bounds on every cut).
% Outputs: studies/gate_density/figures/GD_L_cfp_bounds.png, .csv and
% GD_L_cfp_bounds.mat.

% DEV: the clutter-filter bounds are not in the toolbox. Apply
% studies/sandbox/dev/cfp_error_toolbox.patch first (git apply), and revert
% it after (git apply -R); outputs of the 3 November 2026 run are cached.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");

%% Settings ---------------------------------------------------------------
beam = bmm4d.radar.gate.beam();
table_s = bmm4d.radar.gate.bird_sum_table(beam, 2000, 20000);
s_floor_cm2 = 0.1;
radius_m = 1000;
max_range_m = 80000;
max_elevation_deg = 10;
max_density = 3000;
bin_edges_km = [2:1:10, 12:2:20, 22.5:2.5:50, 55:5:80];
bin_km = (bin_edges_km(1:end - 1) + bin_edges_km(2:end)) / 2;
n_bin = numel(bin_km);
method_name = ["radar output" "estimator" "estimator, CFP bound" "estimator, CFP bound + point clutter" ...
    "estimator, CFP bound + point clutter on batch cuts only"];
n_method = numel(method_name);
n_sweep_max = 20;
effective_radius = 4 / 3 * 6371000;

%% 1. Density per scan, sweep and range bin ---------------------------------------
scans = load(fullfile(bmm4d.util.local_data_dir(), "gate_density", "post2020_scans", "scans.mat")).scans;
n_scan = height(scans);
per_scan = cell(n_scan, 1);
parfor i_s = 1:n_scan
    level2 = bmm4d.radar.read_level2_moments(scans.raw_file(i_s), ["REF" "CFP"], max_elevation_deg);
    level2 = level2(arrayfun(@(s) ~isempty(s.REF), level2));
    % A cut scanned twice (split) keeps its first, surveillance pass, which has CFP
    elevation_key = round([level2.elev] * 2.5) / 2.5;
    is_split = arrayfun(@(k) sum(elevation_key == elevation_key(k)) > 1, 1:numel(level2));
    [~, keep] = unique(elevation_key, "stable");
    level2 = level2(keep);
    is_split = is_split(keep);
    out = struct("density_sum", zeros(n_bin, n_sweep_max, n_method), "gate_count", zeros(n_bin, n_sweep_max, n_method), ...
        "raised", zeros(n_bin, n_sweep_max), "elevation", nan(n_sweep_max, 1), "is_split", false(n_sweep_max, 1));
    for i_e = 1:min(numel(level2), n_sweep_max)
        cut = level2(i_e);
        n_gate = min([size(cut.REF, 1), size(cut.CFP_code, 1), find(cut.REF_range <= max_range_m, 1, "last")]);
        range_m = cut.REF_range(1:n_gate);
        dbz = cut.REF(1:n_gate, :);
        code = double(cut.CFP_code(1:n_gate, :));
        az_step_deg = 360 / size(dbz, 2);
        % CFP: 0 filter not applied, 1 point clutter filter, 2 dual-pol filter, >= 8 removed (code - 8 dB)
        removed_db = code - 8;
        removed_db(code < 8) = nan;
        point_clutter = code == 1;
        output = z_to_refl(idb(dbz)) / 11;
        output(isnan(dbz)) = 0;
        output(output > max_density) = nan;
        detection_dbz = fillmissing(movmedian(prctile(dbz, 1, 2), 9, "omitnan"), "nearest");
        no_echo = isnan(dbz);
        per_gate = cell(1, n_method);
        per_gate{1} = output;
        for i_m = 2:n_method
            switch i_m
                case 2
                    [lo, hi] = bmm4d.radar.gate.bounds(dbz, no_echo, range_m, beam, detection_dbz, s_floor_cm2);
                case 3
                    [lo, hi] = bmm4d.radar.gate.bounds(dbz, no_echo, range_m, beam, detection_dbz, s_floor_cm2, removed_db);
                case 4
                    [lo, hi] = bmm4d.radar.gate.bounds(dbz, no_echo, range_m, beam, detection_dbz, s_floor_cm2, removed_db, point_clutter);
                case 5
                    % Split cuts are filtered only where the radar found clutter: no bound there
                    if is_split(i_e)
                        [lo, hi] = bmm4d.radar.gate.bounds(dbz, no_echo, range_m, beam, detection_dbz, s_floor_cm2);
                    else
                        [lo, hi] = bmm4d.radar.gate.bounds(dbz, no_echo, range_m, beam, detection_dbz, s_floor_cm2, removed_db, point_clutter);
                    end
            end
            estimate = bmm4d.radar.gate.estimate_sweep(lo, hi, range_m, az_step_deg, beam, table_s, radius_m=radius_m);
            per_gate{i_m} = estimate.density;
        end
        bin_of_gate = discretize(range_m / 1000, bin_edges_km);
        rows = ~isnan(bin_of_gate);
        for i_m = 1:n_method
            values = per_gate{i_m}(rows, :);
            valid = ~isnan(values);
            values(~valid) = 0;
            out.density_sum(:, i_e, i_m) = accumarray(bin_of_gate(rows), sum(values, 2), [n_bin 1]);
            out.gate_count(:, i_e, i_m) = accumarray(bin_of_gate(rows), sum(valid, 2), [n_bin 1]);
        end
        raised = per_gate{4}(rows, :) > 2 * per_gate{2}(rows, :);
        out.raised(:, i_e) = accumarray(bin_of_gate(rows), sum(raised, 2), [n_bin 1]);
        out.elevation(i_e) = cut.elev;
        out.is_split(i_e) = is_split(i_e);
    end
    per_scan{i_s} = out;
    fprintf("[cfp bounds] %d/%d %s done\n", i_s, n_scan, scans.station(i_s))
end
% Assemble [n_bin n_sweep n_scan n_method]
parts = cellfun(@(o) o.density_sum, per_scan, UniformOutput=false);
density_sum = permute(cat(4, parts{:}), [1 2 4 3]);
parts = cellfun(@(o) o.gate_count, per_scan, UniformOutput=false);
gate_count = permute(cat(4, parts{:}), [1 2 4 3]);
parts = cellfun(@(o) o.raised, per_scan, UniformOutput=false);
raised_count = cat(3, parts{:});
parts = cellfun(@(o) o.elevation, per_scan, UniformOutput=false);
sweep_elevation = cat(2, parts{:});
parts = cellfun(@(o) o.is_split, per_scan, UniformOutput=false);
sweep_is_split = cat(2, parts{:});
save(fullfile(figure_dir, "GD_L_cfp_bounds.mat"), "density_sum", "gate_count", "raised_count", "sweep_elevation", ...
    "sweep_is_split", "bin_edges_km", "bin_km", "method_name", "scans")

%% 2. Range bias pooled and per group of cuts ----------------------------------------
% Groups: split <= 1.0 deg, split 1.3-2.5 deg, batch 1.8-2.5 deg, batch 3-4.5 deg,
% batch 5-10 deg. Per group, each held at 1 over 30-60 km (all cuts agree far
% from the radar). The pooled fit alone can hide errors that cancel between groups.
n_scan = height(scans);
[i_bin, i_sweep, i_scan] = ndgrid(1:n_bin, 1:n_sweep_max, 1:n_scan);
range_obs_m = bin_km(i_bin(:))' * 1000;
elevation_obs = sweep_elevation(sub2ind(size(sweep_elevation), i_sweep(:), i_scan(:)));
split_obs = sweep_is_split(sub2ind(size(sweep_is_split), i_sweep(:), i_scan(:)));
antenna_obs = scans.antenna_agl(i_scan(:));
group = zeros(size(elevation_obs));
group(split_obs & elevation_obs <= 1.0) = 1;
group(split_obs & elevation_obs > 1.0) = 2;
group(~split_obs & elevation_obs < 2.7) = 3;
group(~split_obs & elevation_obs >= 2.7 & elevation_obs < 4.7) = 4;
group(~split_obs & elevation_obs >= 4.7) = 5;
group_name = ["split ≤ 1.0°" "split 1.3-2.5°" "batch 1.8-2.5°" "batch 3-4.5°" "batch 5-10°"];
low_density = scans.bird_density < 20;
subsets = {true(n_scan, 1), low_density};
subset_name = [sprintf("all %d scans", n_scan) sprintf("%d scans below 20 birds km^{-3}", nnz(low_density))];
methods_shown = [2 4 5];
shown_title = ["estimator (production)" "CFP bound + point clutter, all cuts" "same, batch cuts only"];
pooled = nan(n_bin, n_method, 2);
per_group = nan(n_bin, numel(group_name), n_method, 2);
for i_sub = 1:2
    for i_m = 2:n_method
        mean_density = reshape(density_sum(:, :, :, i_m) ./ gate_count(:, :, :, i_m), [], 1);
        count = reshape(gate_count(:, :, :, i_m), [], 1) .* subsets{i_sub}(i_scan(:));
        fit = gd_fit_range_bias(mean_density, count, range_obs_m, elevation_obs, i_scan(:), antenna_obs, beam, bin_edges_km=bin_edges_km);
        pooled(:, i_m, i_sub) = fit.bias;
        fit = gd_fit_range_bias(mean_density, count, range_obs_m, elevation_obs, i_scan(:), antenna_obs, beam, ...
            bin_edges_km=bin_edges_km, group=group, normalize_groups=true);
        per_group(:, 1:size(fit.bias, 2), i_m, i_sub) = fit.bias;
    end
end
% Absolute level at 30-60 km against the production estimator, by group
far = reshape(bin_km >= 30 & bin_km <= 60, [], 1);
far_level = nan(numel(group_name), n_method);
for i_g = 1:numel(group_name)
    in_group = reshape(group == i_g, n_bin, n_sweep_max, n_scan) & far;
    for i_m = 1:n_method
        far_level(i_g, i_m) = sum(density_sum(:, :, :, i_m) .* in_group, "all") / sum(gate_count(:, :, :, i_m) .* in_group, "all");
    end
end
far_level = far_level ./ far_level(:, 2);
result = table();
for i_sub = 1:2
    for i_m = methods_shown
        result = [result; table(repmat(subset_name(i_sub), n_bin, 1), repmat(method_name(i_m), n_bin, 1), bin_km', ...
            pooled(:, i_m, i_sub), per_group(:, 1, i_m, i_sub), per_group(:, 2, i_m, i_sub), per_group(:, 3, i_m, i_sub), ...
            per_group(:, 4, i_m, i_sub), VariableNames=["scans" "method" "range_km" "pooled" "split_low" "split_high" ...
            "batch_low" "batch_mid"])]; %#ok<AGROW>
    end
end
writetable(result, fullfile(figure_dir, "GD_L_cfp_bounds.csv"))
disp(array2table(far_level, RowNames=group_name, VariableNames=method_name))

%% Plot results -------------------------------------------------------------------
group_colors = [0.15 0.15 0.15; 0.55 0.55 0.55; 0.85 0.37 0.05; 0.95 0.65 0.25];
fig = bmm4d.plot.create_export_figure(3300, 1900, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 2, numel(methods_shown), TileSpacing="tight", Padding="tight");
for i_sub = 1:2
    for k = 1:numel(methods_shown)
        i_m = methods_shown(k);
        ax = nexttile(t);
        hold(ax, "on")
        for i_g = 1:4
            plot(ax, bin_km, per_group(:, i_g, i_m, i_sub), "-", Color=group_colors(i_g, :), LineWidth=1.8, DisplayName=group_name(i_g))
        end
        plot(ax, bin_km, pooled(:, i_m, i_sub), "k--", LineWidth=2.5, DisplayName="all cuts pooled")
        yline(ax, 1, "k:", HandleVisibility="off")
        set(ax, XScale="log", XLim=[2 80], YLim=[0 2.5], FontSize=12, XTick=[2 5 10 20 40 80])
        grid(ax, "on")
        if i_sub == 1
            title(ax, shown_title(k), FontSize=13)
        else
            xlabel(ax, "Distance (km)")
        end
        if k == 1
            ylabel(ax, ["\bf" + subset_name(i_sub) "\rmb(d), each group 1 at 30-60 km"])
        end
        if i_sub == 1 && k == 1
            legend(ax, Location="northeast", FontSize=11)
        end
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_L_cfp_bounds.png"), Resolution=300)
close(fig)
