%% Gate density F: where the radar's clutter filter removes echo power
% Level II CFP (clutter filter power removed) of the sandbox volumes: 28
% radars, frames 2, 4, 6, cuts below 10 deg. CFP codes: 0 filter not
% applied, 1 point clutter filter, >= 8 power removed (code - 8 dB). Per
% class of cut and range bin, the share of gates the filter worked on (any
% code but 0) and the share that lost at least 3 dB (half their power or
% more), among gates with a value or with the filter applied. Classes by cut
% type (split: the elevation is scanned twice, the surveillance pass is
% used; batch: one pass) and elevation. Bins with fewer than MIN_GATES are
% left out. Only gates with the beam centre below MAX_HEIGHT_M above the
% antenna are counted: the layer the application uses, which the higher cuts
% leave within a short range.
% Outputs: studies/gate_density/figures/GD_F_clutter_filter.png, .csv

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");
c = case_settings("sandbox");
manifest = load(fullfile(c.data_dir, "manifest.mat"));

%% Settings ---------------------------------------------------------------
frames = [2 4 6];
max_range_m = 80000;
max_elevation_deg = 10;
max_height_m = 4000;
earth_radius_m = 4 / 3 * 6371000;                  % effective radius, standard refraction
bin_edges_km = [2:1:10, 12:2:20, 22.5:2.5:50, 55:5:80];
min_gates = 5000;
bin_km = (bin_edges_km(1:end - 1) + bin_edges_km(2:end)) / 2;
n_bin = numel(bin_km);
% Classes: cut type, then elevation
elevation_label = ["split 0.5-0.9°" "split 1.3°" "batch 1.8-2.4°" "batch 3.1-5.1°" "batch 6-10°"];
n_elevation = numel(elevation_label);
class_of_cut = @(elev, is_split) (is_split & elev < 1.1) * 1 + (is_split & elev >= 1.1) * 2 + ...
    (~is_split & elev < 2.7) * 3 + (~is_split & elev >= 2.7 & elev < 5.75) * 4 + (~is_split & elev >= 5.75) * 5;

%% Count filtered gates and removed power ---------------------------------------
n_gate = zeros(n_bin, n_elevation);
n_halved = zeros(n_bin, n_elevation);
n_applied = zeros(n_bin, n_elevation);
for i_r = 1:size(manifest.raw_files, 1)
    for i_t = frames
        raw_file = manifest.raw_files(i_r, i_t);
        if ismissing(raw_file) || ~isfile(raw_file)
            continue
        end
        level2 = bmm4d.radar.read_level2_moments(raw_file, ["REF" "CFP"], max_elevation_deg);
        % A split cut repeats its elevation (surveillance pass with CFP, then the Doppler pass)
        all_key = round([level2.elev] * 2.5) / 2.5;
        level2 = level2(arrayfun(@(s) ~isempty(s.CFP) && ~isempty(s.REF), level2));
        elevation_key = round([level2.elev] * 2.5) / 2.5;
        is_split = arrayfun(@(k) sum(all_key == k) > 1, elevation_key);
        [~, keep] = unique(elevation_key, "stable");
        level2 = level2(keep);
        is_split = is_split(keep);
        for i_e = 1:numel(level2)
            i_cut = class_of_cut(level2(i_e).elev, is_split(i_e));
            n_range = min([size(level2(i_e).REF, 1), size(level2(i_e).CFP, 1), find(level2(i_e).REF_range <= max_range_m, 1, "last")]);
            code = double(level2(i_e).CFP_code(1:n_range, :));
            has_value = ~isnan(level2(i_e).REF(1:n_range, :));
            filtered = code >= 8;
            halved = code >= 11;                     % at least 3 dB removed
            range_m = level2(i_e).REF_range(1:n_range);
            height_m = range_m * sind(level2(i_e).elev) + range_m .^ 2 / (2 * earth_radius_m);
            bin_of_gate = discretize(range_m / 1000, bin_edges_km);
            bin_of_gate(height_m > max_height_m) = nan;
            for k = 1:n_bin
                rows = bin_of_gate == k;
                if ~any(rows)
                    continue
                end
                considered = has_value(rows, :) | filtered(rows, :);
                n_gate(k, i_cut) = n_gate(k, i_cut) + nnz(considered);
                n_halved(k, i_cut) = n_halved(k, i_cut) + nnz(halved(rows, :));
                n_applied(k, i_cut) = n_applied(k, i_cut) + nnz(filtered(rows, :) | code(rows, :) == 1);
            end
        end
    end
    fprintf("[clutter filter] radar %d done\n", i_r)
end

%% Share of gates that lost half their power or more -------------------------------
share_filtered = n_halved ./ n_gate;
share_filtered(n_gate < min_gates) = nan;
share_applied = n_applied ./ n_gate;
share_applied(n_gate < min_gates) = nan;
summary = table(bin_km', VariableNames="range_km");
for i_c = 1:n_elevation
    summary.(matlab.lang.makeValidName("lost3dB_" + elevation_label(i_c))) = round(share_filtered(:, i_c), 2);
    summary.(matlab.lang.makeValidName("applied_" + elevation_label(i_c))) = round(share_applied(:, i_c), 2);
end
disp(summary)
writetable(summary, fullfile(figure_dir, "GD_F_clutter_filter.csv"))

%% Plot: where the filter works, and where it takes half the power -------------------------
% Split cuts in greys, batch cuts in oranges (as the cut-level figures)
cut_colors = [0.15 0.15 0.15; 0.55 0.55 0.55; 0.8 0.3 0.0; 0.92 0.55 0.15; 0.98 0.78 0.45];
fig = bmm4d.plot.create_export_figure(3000, 1100, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
shares = {share_applied, share_filtered};
titles = ["Gates the filter worked on" "Gates that lost ≥ 3 dB"];
for k = 1:2
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    for i_c = 1:n_elevation
        plot(ax, bin_km, 100 * shares{k}(:, i_c), "-o", Color=cut_colors(i_c, :), LineWidth=2, MarkerSize=4, ...
            MarkerFaceColor=cut_colors(i_c, :), DisplayName=elevation_label(i_c))
    end
    set(ax, XScale="log", FontSize=12); xlim(ax, [2 80]); ylim(ax, [0 100]); xticks(ax, [2 5 10 20 40 80])
    xlabel(ax, "Range (km)"); title(ax, titles(k), FontWeight="normal")
    if k == 1
        ylabel(ax, "Share of gates (%)")
        legend(ax, Location="southwest", FontSize=11)
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "GD_F_clutter_filter.png"), Resolution=150)
close(fig)
