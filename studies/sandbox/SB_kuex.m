%% Sandbox: KUEX on 21-22 October 2026, the lowest 1000 m within 30 km
% On the observation sections (SB_obs_sections) KUEX on this night is a
% clear case of sweeps disagreeing near the radar: below about 600 m the
% upper sweeps read a few birds km^-3 and the lowest ones 20-40. What is
% in the data there?
%   1. Maps of each sweep within 30 km: radar output (REF), power removed
%      by the clutter filter (CFP), radial velocity and the cleaned
%      density that enters the inversion (gate estimator, full gates).
%   2. Range-height section of the azimuth means of the same quantities,
%      with each sweep's beam (half-power) as a band.
%   3. Beam consistency: a vertical profile from the narrow beams only
%      (sweeps >= 1.8 deg within 25 km: beam <= 400 m) passed through each
%      lower sweep's beam (0.925 deg, as the forward model) at its range,
%      against what that sweep measured. If the near sweeps are right and
%      the low layer is nearly empty, the low sweeps should read what the
%      smeared profile gives; if they read more, they carry echo the narrow
%      beams do not see (clutter, wider beam than modelled, refraction).
%   Times: 02, 04 and 06 UTC on 22 October 2026.
% Outputs: figures/kuex_maps.png, kuex_section.png, kuex_beam.png.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
c = case_settings("night_20261021");
cache_dir = bmm4d.radar.volume_cache_dir(c);
files = dir(fullfile(cache_dir, "KUEX20261022_*.mat"));
time = datetime(extractBetween(string({files.name}), 5, 19), InputFormat="yyyyMMdd_HHmmss");
target = datetime(2026, 10, 22, [2 4 6], 0, 0);
[~, pick] = min(abs(time(:) - target), [], 1);
files = files(pick);
read_options = [c.read_options {"maxRange", 40000, "data_res", [0 0]}];
hpbw = 0.925;                            % forward model half-power beam width (deg)

%% Read raw moments and cleaned density -------------------------------------------
volumes = struct("time", {}, "sweeps", {});
for i_f = 1:numel(files)
    raw_file = load(fullfile(files(i_f).folder, files(i_f).name)).volume.source_file;
    moments = bmm4d.radar.read_level2_moments(raw_file, ["REF" "VEL" "CFP"], 7);
    cleaned = bmm4d.radar.read_rd(raw_file, read_options{:});
    cleaned = cleaned{1};
    % One pass per elevation: the surveillance pass (REF and CFP) of split
    % cuts, with the velocity of its Doppler pass
    sweeps = struct("elev", {}, "az", {}, "range_m", {}, "dbz", {}, "cfp", {}, "vel", {}, "density", {});
    for e = unique(round([moments.elev] * 10) / 10)
        passes = moments(abs([moments.elev] - e) < 0.1 & arrayfun(@(m) ~isempty(m.REF), moments));
        [~, k] = max(arrayfun(@(m) numel(m.REF_range), passes));
        main = passes(k);
        with_vel = passes(arrayfun(@(m) ~isempty(m.VEL), passes));
        n_gate = find(main.REF_range <= 40000, 1, "last");
        vel = nan(n_gate, numel(main.az));
        if ~isempty(with_vel)
            n_vel = min(n_gate, numel(with_vel(1).VEL_range));
            vel(1:n_vel, :) = with_vel(1).VEL(1:n_vel, :);
        end
        cfp = nan(n_gate, numel(main.az));
        if ~isempty(main.CFP), cfp = main.CFP(1:n_gate, :); end
        dbz = main.REF(1:n_gate, :);
        dbz(main.REF_code(1:n_gate, :) == 0) = -Inf;
        % Cleaned density on the same grid (nearest azimuth)
        [~, i_clean] = min(abs([cleaned.dens.sweeps.elev] - main.elev));
        clean_sweep = cleaned.dens.sweeps(i_clean);
        [~, az_id] = min(abs(mod(main.az(:) - clean_sweep.az(:)' + 180, 360) - 180), [], 2);
        range_id = interp1(clean_sweep.range(:), 1:numel(clean_sweep.range), main.REF_range(1:n_gate), "nearest", NaN);
        density = nan(n_gate, numel(main.az));
        ok = ~isnan(range_id);
        density(ok, :) = clean_sweep.data(range_id(ok), az_id);
        sweeps(end + 1) = struct("elev", main.elev, "az", main.az, "range_m", main.REF_range(1:n_gate), ... %#ok<AGROW>
            "dbz", dbz, "cfp", cfp, "vel", vel, "density", density);
    end
    volumes(i_f) = struct("time", time(pick(i_f)), "sweeps", sweeps);
end

%% Plot: maps within 30 km (04 UTC) ----------------------------------------------------
volume = volumes(2);
shown = find(ismember(round([volume.sweeps.elev] * 10) / 10, [0.5 0.9 1.3 1.8 3.1 6.4]));
quantities = ["Radar output REF (dBZ)" "Clutter filter power removed (dB)" "Radial velocity (m s^{-1})" ...
    "Cleaned density (birds km^{-3}, log10)"];
fig = bmm4d.plot.create_export_figure(3200, 800 * numel(shown), DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, numel(shown), numel(quantities), TileSpacing="tight", Padding="tight");
title(t, "KUEX " + string(volume.time, "dd MMM yyyy HH:mm") + " UTC, within 30 km; circles every 10 km", FontSize=16, Color="k")
for i_s = shown
    sweep = volume.sweeps(i_s);
    [az_grid, range_grid] = meshgrid(sweep.az, sweep.range_m / 1000);
    x = range_grid .* sind(az_grid);
    y = range_grid .* cosd(az_grid);
    values = {max(sweep.dbz, -20), sweep.cfp, sweep.vel, log10(max(sweep.density, 0.1))};
    limits = {[-20 30], [0 30], [-20 20], [-1 2.5]};
    maps = {parula, flipud(gray(256)), redblue_map(), parula};
    for i_q = 1:numel(quantities)
        ax = nexttile(t);
        surface(ax, x, y, zeros(size(x)), values{i_q}, EdgeColor="none")
        view(ax, 2); axis(ax, "equal"); xlim(ax, [-30 30]); ylim(ax, [-30 30])
        colormap(ax, maps{i_q}); clim(ax, limits{i_q}); set(ax, Color=[0.9 0.9 0.9], FontSize=10)
        hold(ax, "on")
        for r = 10:10:30
            plot(ax, r * sind(0:360), r * cosd(0:360), "k:", LineWidth=0.5)
        end
        if i_s == shown(1)
            title(ax, quantities(i_q), FontWeight="bold", FontSize=12)
            cb = colorbar(ax, "eastoutside"); cb.FontSize = 9;
        end
        if i_q == 1
            ylabel(ax, sprintf("%.1f°", sweep.elev), FontWeight="bold", FontSize=13)
        end
    end
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "kuex_maps.png"), Resolution=150)
close(fig)

%% Range-height section: azimuth means -----------------------------------------------
% Mean over azimuth per gate range: radar output (no echo = 0), output plus
% the power the clutter filter removed (upper bound: it also adds ground
% clutter), and the cleaned density; one band per sweep (half-power beam)
range_edges = (2000:500:30000)';
fig = bmm4d.plot.create_export_figure(3600, 1100 * numel(volumes), DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, numel(volumes), 3, TileSpacing="tight", Padding="tight");
title(t, "KUEX 22 Oct 2026: azimuth means within 30 km, each sweep drawn over its half-power beam", FontSize=16, Color="k")
panel_names = ["Radar output" "Output + clutter-filter power" "Cleaned density (inversion input)"];
for i_v = 1:numel(volumes)
    for i_q = 1:3
        ax = nexttile(t); hold(ax, "on"); box(ax, "on")
        for sweep = volumes(i_v).sweeps
            z = 10 .^ (sweep.dbz / 10);
            if i_q == 2
                z = z .* 10 .^ (fillmissing(sweep.cfp, "constant", 0) / 10);
            end
            density = z_to_refl(z) / 11;
            if i_q == 3, density = sweep.density; end
            range_bin = discretize(sweep.range_m, range_edges);
            ok = ~isnan(range_bin);
            mean_density = accumarray(range_bin(ok), mean(density(ok, :), 2, "omitnan"), [numel(range_edges) - 1 1], @mean, NaN);
            centre = (range_edges(1:end - 1) + range_edges(2:end)) / 2;
            [~, height_m] = slant2ground(centre, sweep.elev);
            half = centre * tand(hpbw / 2);
            % One patch per range bin, coloured by log10 density
            for i_b = find(isfinite(mean_density))'
                patch(ax, centre(i_b) / 1000 + [-0.25 0.25 0.25 -0.25], height_m(i_b) + half(i_b) * [-1 -1 1 1], ...
                    log10(max(mean_density(i_b), 0.1)), EdgeColor="none", FaceAlpha=0.8)
            end
            plot(ax, centre / 1000, height_m, "k-", LineWidth=0.3)
        end
        colormap(ax, parula); clim(ax, [-1 2.5]); xlim(ax, [0 30]); ylim(ax, [0 1500]); set(ax, FontSize=12)
        title(ax, panel_names(i_q) + ", " + string(volumes(i_v).time, "HH:mm") + " UTC", FontWeight="normal")
        if i_q == 1, ylabel(ax, "Height above antenna (m)"), end
        if i_v == numel(volumes), xlabel(ax, "Range (km)"), end
    end
end
cb = colorbar(ax); cb.Layout.Tile = "east"; cb.Label.String = "log_{10} birds km^{-3}"; cb.FontSize = 12;
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "kuex_section.png"), Resolution=150)
close(fig)

%% Beam consistency ------------------------------------------------------------------
% Profile from the narrow beams (>= 1.8 deg, <= 25 km), 50 m bins, cleaned
% density; then each sweep's prediction = profile smeared by its Gaussian
% beam at each range, against its measured azimuth mean
fig = bmm4d.plot.create_export_figure(3600, 1100, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, numel(volumes), TileSpacing="tight", Padding="tight");
title(t, "KUEX: what each sweep measures (solid) and what the narrow-beam profile predicts through its beam (dashed)", ...
    FontSize=16, Color="k")
profile_edges = (0:50:2500)';
profile_centre = (profile_edges(1:end - 1) + profile_edges(2:end)) / 2;
colors = parula(8);
for i_v = 1:numel(volumes)
    heights = []; values = [];
    for sweep = volumes(i_v).sweeps
        if sweep.elev < 1.7, continue, end
        near = sweep.range_m <= 25000 & sweep.range_m >= 3000;
        [~, height_m] = slant2ground(sweep.range_m(near), sweep.elev);
        heights = [heights; height_m]; %#ok<AGROW>
        values = [values; mean(sweep.density(near, :), 2, "omitnan")]; %#ok<AGROW>
    end
    bin = discretize(heights, profile_edges);
    ok = ~isnan(bin) & ~isnan(values);
    narrow_profile = accumarray(bin(ok), values(ok), [numel(profile_centre) 1], @mean, NaN);
    narrow_profile = fillmissing(narrow_profile, "linear", EndValues="nearest");
    ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
    low = find([volumes(i_v).sweeps.elev] < 1.7);
    for i_l = 1:numel(low)
        sweep = volumes(i_v).sweeps(low(i_l));
        range_m = (3000:500:40000)';
        [~, centre_m] = slant2ground(range_m, sweep.elev);
        sd = range_m * deg2rad(hpbw) / (2 * sqrt(2 * log(2)));
        weight = exp(-(profile_centre' - centre_m) .^ 2 ./ (2 * sd .^ 2));
        predicted = (weight * narrow_profile) ./ sum(weight, 2);
        range_bin = discretize(sweep.range_m, [range_m - 250; range_m(end) + 250]);
        ok = ~isnan(range_bin);
        measured = accumarray(range_bin(ok), mean(sweep.density(ok, :), 2, "omitnan"), [numel(range_m) 1], @mean, NaN);
        color = colors(2 * i_l - 1, :);
        plot(ax, range_m / 1000, measured, "-", Color=color, LineWidth=2, DisplayName=sprintf("%.1f° measured", sweep.elev))
        plot(ax, range_m / 1000, predicted, "--", Color=color, LineWidth=2, DisplayName=sprintf("%.1f° from narrow beams", sweep.elev))
    end
    set(ax, YScale="log", FontSize=12); xlim(ax, [0 40]); ylim(ax, [0.3 300])
    xlabel(ax, "Range (km)"); ylabel(ax, "Azimuth-mean density (birds km^{-3})")
    title(ax, string(volumes(i_v).time, "HH:mm") + " UTC", FontWeight="normal")
    legend(ax, Location="southeast", FontSize=10)
    % Inset: the narrow-beam profile
    inset = axes(fig, Position=ax.Position .* [1 1 0.3 0.4] + [ax.Position(3) * 0.08 ax.Position(4) * 0.55 0 0]);
    plot(inset, narrow_profile, profile_centre, "k-", LineWidth=1.5)
    set(inset, XScale="log", FontSize=9); ylim(inset, [0 2000]); xlim(inset, [0.3 300])
    title(inset, "narrow-beam profile", FontSize=9, FontWeight="normal")
end
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "kuex_beam.png"), Resolution=150)
close(fig)

%% Local functions ----------------------------------------------------------
function map = redblue_map()
% Diverging blue-white-red for radial velocity
n = 128;
map = [linspace(0.2, 1, n)' linspace(0.3, 1, n)' ones(n, 1); ones(n, 1) linspace(1, 0.25, n)' linspace(1, 0.2, n)'];
end
