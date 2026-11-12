%% Sandbox: the observations on a range-height section, per radar and night
% The cleaned, thinned gates that enter the inversion (gate estimator
% density, volume cache of the sandbox cleaning settings), shown as one dot
% per sweep and 0.5 km range bin at the beam centre height, coloured by
% the median density over all azimuths (and all volumes of the night where
% the night was run in full). Where sweeps cross the same height, the dots
% should have the same colour.
% Nights: the two full nights (night_20261021, night_20210512; 16-23
% volumes per radar) and the nine dry nights of SB_night_screen with one
% volume per radar about 3 h after sunset (sandbox_nights/raw, cleaned here
% with the production settings). Five of those volumes were taken the
% evening before the screened night (getSunriseSunset bug, MIGRATION.md):
% nights are labelled by the volume's own date.
% Score per radar and night: median over sweep pairs of log2(upper sweep /
% lower sweep) at the same height, upper sweep within 3-20 km (0 = the
% sweeps agree, < 0 = the near-radar deficit). The two best and two worst
% radars (mean over nights) are plotted.
% Outputs: figures/obs_sections.png, obs_sections_score.csv,
% obs_sections_density.png and .csv (score against the density aloft).

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
range_edges_m = 3000:500:60000;
range_centre_m = (range_edges_m(1:end - 1) + range_edges_m(2:end)) / 2;
c = case_settings("sandbox");
cache_dir = bmm4d.radar.volume_cache_dir(c);

%% Clean the single volumes of the screened nights -------------------------------
raw_list = dir("/Users/rafnuss/Data/BMM-within/sandbox_nights/raw/K*_V06");
raw_files = string(fullfile({raw_list.folder}, {raw_list.name}))';
out_files = fullfile(cache_dir, string({raw_list.name})' + ".mat");
todo = find(~isfile(out_files));
fprintf("%d volumes to clean\n", numel(todo))
read_options = c.read_options; vp = c.vp; filter = c.filter;
parfor i_t = 1:numel(todo)
    bmm4d.radar.clean_volume(raw_files(todo(i_t)), read_options, vp, filter, out_files(todo(i_t)));
end

%% Volumes of each radar and night ------------------------------------------------
volumes = table(out_files(isfile(out_files)), VariableNames="file");
for night = ["night_20261021" "night_20210512"]
    c_night = case_settings(night);
    files = dir(fullfile(cache_dir, "K*.mat"));
    names = string({files.name})';
    time = datetime(extractBetween(names, 5, 19), InputFormat="yyyyMMdd_HHmmss");
    keep = time >= c_night.times(1) - minutes(15) & time <= c_night.times(end) + minutes(15);
    volumes = [volumes; table(fullfile(cache_dir, names(keep)), VariableNames="file")]; %#ok<AGROW>
end
volumes = unique(volumes);
[~, names] = fileparts(volumes.file);
volumes.radar = extractBefore(names, 5);
volumes.time = datetime(extractBetween(names, 5, 19), InputFormat="yyyyMMdd_HHmmss");
volumes.night = string(dateshift(volumes.time - hours(12), "start", "day"), "yyyy-MM-dd");
nights = unique(volumes.night);
radars = unique(volumes.radar);
fprintf("%d volumes, %d nights, %d radars\n", height(volumes), numel(nights), numel(radars))

%% Median density per sweep and range bin, and the score ---------------------------
sections = cell(numel(radars), numel(nights));
score = nan(numel(radars), numel(nights));
for i_n = 1:numel(nights)
    for i_r = 1:numel(radars)
        files = volumes.file(volumes.radar == radars(i_r) & volumes.night == nights(i_n));
        if isempty(files), continue, end
        % Pool the kept gates of every volume by sweep (elevation to 0.1 deg) and range bin
        pooled = zeros(0, 3);
        for file = files'
            volume = load(file).volume;
            for sweep = volume.dens.sweeps
                range_bin = discretize(sweep.range(:), range_edges_m);
                kept = ~isnan(sweep.data) & ~isnan(range_bin);
                [row, ~] = find(kept);
                pooled = [pooled; repmat(round(sweep.elev * 10) / 10, numel(row), 1) range_bin(row) sweep.data(kept)]; %#ok<AGROW>
            end
        end
        [sweep_id, elevations] = findgroups(pooled(:, 1));
        section = table();
        for i_e = 1:numel(elevations)
            in_sweep = sweep_id == i_e;
            median_density = accumarray(pooled(in_sweep, 2), pooled(in_sweep, 3), [numel(range_centre_m) 1], @median, NaN);
            n_gate = accumarray(pooled(in_sweep, 2), 1, [numel(range_centre_m) 1]);
            [~, height_m] = slant2ground(range_centre_m', elevations(i_e));
            ok = n_gate >= 10;
            section = [section; table(repmat(elevations(i_e), nnz(ok), 1), range_centre_m(ok)' / 1000, height_m(ok), ...
                median_density(ok), VariableNames=["elev" "range_km" "height_m" "density"])]; %#ok<AGROW>
        end
        sections{i_r, i_n} = section;
        % Score: each upper-sweep dot within 3-20 km against the next lower
        % sweep interpolated (log density) to the same height
        log_ratio = [];
        for i_e = 2:numel(elevations)
            upper = section(section.elev == elevations(i_e) & section.range_km <= 20 & section.density > 0.5, :);
            lower = section(section.elev == elevations(i_e - 1) & section.density > 0.5, :);
            if height(upper) == 0 || height(lower) < 2, continue, end
            lower_at_height = interp1(lower.height_m, log2(lower.density), upper.height_m);
            log_ratio = [log_ratio; log2(upper.density) - lower_at_height]; %#ok<AGROW>
        end
        if nnz(isfinite(log_ratio)) >= 5
            score(i_r, i_n) = median(log_ratio, "omitnan");
        end
    end
end
score_table = array2table(round(score, 2), VariableNames=nights, RowNames=radars);
score_table.mean = round(mean(score, 2, "omitnan"), 2);
score_table = sortrows(score_table, "mean");
disp(score_table)
writetable(score_table, fullfile(figure_dir, "obs_sections_score.csv"), WriteRowNames=true)

%% Plot results -----------------------------------------------------------
% Two worst and two best radars (columns), every night (rows). Colour:
% log10 density, shared limits; below 0.5 birds km^-3 at the bottom
shown = [string(score_table.Properties.RowNames(1:2)); string(score_table.Properties.RowNames(end - 1:end))];
color_limits = [-0.3 2.5];
fig = bmm4d.plot.create_export_figure(3200, 6000, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, numel(nights), numel(shown), TileSpacing="tight", Padding="tight");
title(t, "Observations entering the inversion: two worst (left) and two best (right) radars, median over azimuth", ...
    FontSize=17, Color="k")
xlabel(t, "Range (km)", FontSize=15)
ylabel(t, "Height above antenna (m)", FontSize=15)
for i_n = 1:numel(nights)
    for i_s = 1:numel(shown)
        i_r = find(radars == shown(i_s));
        ax = nexttile(t); hold(ax, "on"); box(ax, "on")
        section = sections{i_r, i_n};
        if isempty(section), axis(ax, "off"), continue, end
        for e = unique(section.elev)'
            in_sweep = section.elev == e;
            plot(ax, section.range_km(in_sweep), section.height_m(in_sweep), "-", Color=[0.85 0.85 0.85], LineWidth=0.5)
        end
        scatter(ax, section.range_km, section.height_m, 14, log10(max(section.density, 0.5)), "filled")
        colormap(ax, parula); clim(ax, color_limits)
        xlim(ax, [0 60]); ylim(ax, [0 2500]); set(ax, FontSize=11, Color="w")
        % Radar and mean score as column headers, night as row label, the panel's score inside
        n_volume = nnz(volumes.radar == shown(i_s) & volumes.night == nights(i_n));
        if i_n == 1
            title(ax, sprintf("%s, mean %.2f",shown(i_s), score_table{shown(i_s), "mean"}), FontWeight="bold", FontSize=14)
        end
        if i_s == 1
            ylabel(ax, sprintf("%s\n%d vol.", nights(i_n), n_volume), FontWeight="bold", FontSize=12)
        end
        text(ax, 58, 2400, sprintf("%.2f", score(i_r, i_n)), HorizontalAlignment="right", VerticalAlignment="top", ...
            FontSize=12, FontWeight="bold", BackgroundColor="w")
    end
end
cb = colorbar(ax);
cb.Layout.Tile = "east";
cb.Label.String = "Median density (birds km^{-3})";
cb.Ticks = log10([0.5 1 3 10 30 100 300]);
cb.TickLabels = ["≤0.5" "1" "3" "10" "30" "100" "300"];
cb.FontSize = 13;
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "obs_sections.png"), Resolution=150)
close(fig)

%% Score against the density aloft ------------------------------------------------
% Median of the section dots at 500-1500 m per radar-night
density_aloft = nan(size(score));
for i_s = 1:numel(sections)
    if isempty(sections{i_s}), continue, end
    aloft = sections{i_s}.height_m >= 500 & sections{i_s}.height_m < 1500;
    density_aloft(i_s) = median(sections{i_s}.density(aloft), "omitnan");
end
ok = isfinite(score) & density_aloft > 0;
writetable(table(density_aloft(ok), score(ok), VariableNames=["density_aloft" "score"]), ...
    fullfile(figure_dir, "obs_sections_density.csv"))
fprintf("Spearman correlation of score and log density aloft: %.2f (%d radar-nights)\n", ...
    corr(score(ok), log(density_aloft(ok)), Type="Spearman"), nnz(ok))
night_index = repmat(1:numel(nights), numel(radars), 1);
fig = bmm4d.plot.create_export_figure(2000, 1300, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
ax = axes(fig, Position=[0.1 0.12 0.72 0.8]); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
scatter(ax, density_aloft(ok), score(ok), 40, night_index(ok), "filled", MarkerFaceAlpha=0.8)
colormap(ax, copper(numel(nights)));
cb = colorbar(ax); cb.Ticks = 1:numel(nights); cb.TickLabels = nights; cb.FontSize = 10;
yline(ax, 0, "k:")
set(ax, XScale="log", FontSize=13); xlim(ax, [0.1 1000]); ylim(ax, [-2 0.6])
xlabel(ax, "Median density at 500-1500 m (birds km^{-3})")
ylabel(ax, "Score: log_2 upper / lower sweep")
title(ax, sprintf("28 radars x 11 nights (Spearman %.2f)", ...
    corr(score(ok), log(density_aloft(ok)), Type="Spearman")), FontWeight="normal", FontSize=15)
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "obs_sections_density.png"), Resolution=150)
close(fig)
