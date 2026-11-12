%% Sandbox: the clutter filter's removed power as an observation error
% The radar's clutter filter (Level II CFP) removes power from bird echo as
% well as from ground clutter, in shares that change with range, cut and
% density (gate-density report Part III, gate_density GD_L). Moving
% the estimate up by the removed power lets clutter in. Here the estimate
% stays, and the filter only lowers the weight of what it touched: SD_CFP
% (bmm4d.radar.filter.estimate_density, read_rd option 'cfp') is the
% estimate with each filtered gate taken as at most the output plus the
% removed power, minus the production estimate. It is added in quadrature
% to the observation error (bmm4d.solve.observation_weights).
% 1. Clean the 168 sandbox volumes again with CFP (prod/volumes/sandbox_cfp)
%    and write copies whose data is SD_CFP;
% 2. observation rows of SD_CFP (prod/sandbox_cfp), the same rows as the
%    density;
% 3. solves: production weights (estimator sd at the prior, floor 1.5) and
%    the same plus SD_CFP;
% 4. sections through KDMX, KLSX and KDVN and column maps at every time.
% Outputs: figures/cfp_error_weights.csv, cfp_error_sections.gif and
% cfp_error_sections_<HHmm>.png.

% DEV: the clutter-filter bounds are not in the toolbox. Apply
% studies/sandbox/dev/cfp_error_toolbox.patch first (git apply), and revert
% it after (git apply -R); outputs of the 3 November 2026 run are cached.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "sandbox"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
S = sandbox_load("sandbox");
c = case_settings("sandbox");
manifest = load(fullfile(c.data_dir, "manifest.mat"));
cfp_dir = fullfile(bmm4d.util.local_data_dir(), "prod", "volumes", "sandbox_cfp");
density_dir = fullfile(cfp_dir, "density");
error_dir = fullfile(cfp_dir, "sd_cfp");

%% 1. Clean again with CFP ---------------------------------------------------------
[~, stem] = fileparts(manifest.volume_files);
has = manifest.volume_files ~= "";
density_files = strings(size(manifest.volume_files));
density_files(has) = fullfile(density_dir, stem(has) + ".mat");
error_files = strings(size(manifest.volume_files));
error_files(has) = fullfile(error_dir, stem(has) + ".mat");
if ~isfolder(density_dir), mkdir(density_dir); end
if ~isfolder(error_dir), mkdir(error_dir); end
% Volumes cleaned already are kept (resumable)
todo = find(has & ~isfile(density_files));
raw_files = manifest.raw_files(todo);
out_files = density_files(todo);
read_options = [c.read_options {"cfp", true}];
vp = c.vp;
filter = c.filter;
filter.cfp_error = true;
status = strings(numel(todo), 1);
if ~isempty(todo)
    % Four workers: each holds the bird-sum table and a full volume, and 12 ran out of memory
    pool = gcp("nocreate");
    if isempty(pool) || pool.NumWorkers > 4
        delete(pool);
        parpool("Processes", 4);
    end
    parfor k = 1:numel(todo)
        status(k) = bmm4d.radar.clean_volume(raw_files(k), read_options, vp, filter, out_files(k));
    end
end
disp(unique(status))
% Copy with SD_CFP as the data (the rows of step 2), and check the density
% against the production cache
max_difference = 0;
for k = find(has & ~isfile(error_files))'
    volume = load(density_files(k)).volume;
    cached = load(manifest.volume_files(k)).volume;
    for i_e = 1:numel(volume.dens.sweeps)
        max_difference = max(max_difference, max(abs(volume.dens.sweeps(i_e).data - cached.dens.sweeps(i_e).data), [], "all", "omitnan"));
        volume.dens.sweeps(i_e).data = volume.dens.sweeps(i_e).sd_cfp;
    end
    volume.dens.sweeps = rmfield(volume.dens.sweeps, intersect(["sd" "sd_cfp"], string(fieldnames(volume.dens.sweeps))));
    save(error_files(k), "volume")
end
fprintf("Largest difference to the cached density: %.3g birds km^-3\n", max_difference)

%% 2. Observation rows of SD_CFP ---------------------------------------------------------
c_cfp = c;
c_cfp.name = "sandbox_cfp";
c_cfp.data_dir = fullfile(bmm4d.util.local_data_dir(), "prod", c_cfp.name);
if ~isfile(fullfile(c_cfp.data_dir, "case.mat"))
    if ~isfolder(c_cfp.data_dir), mkdir(c_cfp.data_dir); end
    volume_files = error_files; raw_files = manifest.raw_files; radars = manifest.radars; times = manifest.times;
    save(fullfile(c_cfp.data_dir, "manifest.mat"), "raw_files", "radars", "times", "volume_files")
    copyfile(fullfile(c.data_dir, "static_clutter.mat"), c_cfp.data_dir)
    bmm4d.build_observations(c_cfp);
end
loaded = load(fullfile(c_cfp.data_dir, "case.mat"), "observations");
filter_error = loaded.observations.values;
clear loaded
assert(isequal(size(filter_error), size(S.observations.values)), "The SD_CFP rows differ from the sandbox rows.")
% Only where there is a density; the eps of zero values goes
filter_error = filter_error .* (S.observations.values ~= 0);
filter_error(filter_error < 1e-6) = 0;
fprintf("Observations with a filter error > 0: %.1f%%\n", 100 * nnz(filter_error) / nnz(S.observations.values))

%% 3. Solves ---------------------------------------------------------------------------
prediction = S.observations.F(:, S.active_id) * double(sandbox_prior(S));
floor_birds = c.inversion.observation_error_floor;
weights_sd = bmm4d.solve.observation_weights(S.observations, prediction, floor_birds);
weights_cfp = bmm4d.solve.observation_weights(S.observations, prediction, floor_birds, filter_error);
runs = {sandbox_run(S, struct("name", "production weights (sd, floor 1.5)", "video", false, "row_weight", @(r, e, z) weights_sd)), ...
    sandbox_run(S, struct("name", "production weights + clutter filter error", "video", false, "row_weight", @(r, e, z) weights_cfp))};
% Weight relative to production by cut type and range (median over observations)
[row, col] = find(S.observations.values);
linear_id = sub2ind(size(weights_sd), row, col);
range_km = double(S.observations.row_range_m(row)) / 1000;
elevation = double(S.observations.row_elev_deg(row));
cut_type = repmat("batch ≥ 1.8°", size(row));
cut_type(elevation < 0.7) = "split 0.5°";
cut_type(elevation >= 0.7 & elevation < 1.6) = "split 0.9-1.3°";
range_edges_km = [0 10 20 30 40 60 80 120 250];
range_bin = discretize(range_km, range_edges_km);
weight_ratio = weights_cfp(linear_id) ./ weights_sd(linear_id);
summary = groupsummary(table(cut_type, range_bin, weight_ratio), ["cut_type" "range_bin"], "median", "weight_ratio");
summary.range_km = string(range_edges_km(summary.range_bin))' + "-" + string(range_edges_km(summary.range_bin + 1))';
writetable(summary(:, ["cut_type" "range_km" "GroupCount" "median_weight_ratio"]), fullfile(figure_dir, "cfp_error_weights.csv"))
disp(unstack(summary(:, ["cut_type" "range_km" "median_weight_ratio"]), "median_weight_ratio", "range_km"))
for k = 1:2
    fprintf("%s: %.1f s, %d iterations, birds %.3f of the first\n", runs{k}.name, runs{k}.seconds, runs{k}.iterations, ...
        sum(runs{k}.field, "all") / sum(runs{1}.field, "all"))
end

%% 4. Sections and column maps -----------------------------------------------------------------
g = S.g;
radar_selected = S.radar_selected;
section_radars = ["KDMX" "KLSX" "KDVN"];
row_names = ["production weights" "+ clutter filter error" "log_2 ratio"];
effective_radius = 6371200 * 4 / 3;
signed_map = interp1([0 0.5 1], [0.13 0.4 0.75; 1 1 1; 0.8 0.2 0.15], linspace(0, 1, 256));
[ix, iy, ~] = ind2sub(g.sz, S.active_id);
column_id = sub2ind(g.sz(1:2), ix, iy);
air = bmm4d.grid.air_fraction(g, S.active_id);
to_column = @(x) accumarray(column_id, double(x) .* air, [prod(g.sz(1:2)) 1]) * g.dz / 1000;
active_column = to_column(ones(numel(S.active_id), 1)) > 0;
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
n_time = size(runs{1}.field, 2);
gif_file = fullfile(figure_dir, "cfp_error_sections.gif");
fig = bmm4d.plot.create_export_figure(3600, 2400, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
for i_t = 1:n_time
    clf(fig)
    t = tiledlayout(fig, 3, 1 + numel(section_radars), TileSpacing="tight", Padding="tight");
    title(t, "Sandbox " + string(c.times(i_t), "dd MMM yyyy HH:mm") + " UTC (dashed: 0.5° beam centre)", FontSize=18, FontWeight="bold")
    volumes = cell(1, 2);
    for k = 1:2
        volumes{k} = nan(g.sz, "single");
        volumes{k}(S.active_id) = runs{k}.field(:, i_t);
    end
    for i_row = 1:3
        % Column from above
        ax = nexttile(t);
        column = nan(g.sz(1:2));
        if i_row == 3
            ratio = log2(max(to_column(runs{2}.field(:, i_t)), 1) ./ max(to_column(runs{1}.field(:, i_t)), 1));
            column(active_column) = ratio(active_column);
            bmm4d.plot.plot_night_map(ax, g, background, column, [-0.5 0.5], states, radar_selected, c.figure_bounds, ones(g.sz(1:2)));
            colormap(ax, signed_map)
        else
            values = to_column(runs{i_row}.field(:, i_t));
            column(active_column) = log10(max(values(active_column), 1));
            bmm4d.plot.plot_night_map(ax, g, background, column, log10([5 300]), states, radar_selected, c.figure_bounds, ones(g.sz(1:2)));
            colormap(ax, parula)
        end
        ylabel(ax, row_names(i_row), FontWeight="bold", FontSize=14)
        % West-east sections through each radar
        for i_k = 1:numel(section_radars)
            k = find(radar_selected.name == section_radars(i_k), 1);
            [~, y_id] = min(abs(g.lat - radar_selected.lat(k)));
            if i_row == 3
                field = log2(max(volumes{2}, 1) ./ max(volumes{1}, 1));
            else
                field = volumes{i_row};
            end
            values = squeeze(field(:, y_id, :))';
            distance_km = (g.lon - radar_selected.lon(k)) * 111.32 * cosd(radar_selected.lat(k));
            ax = nexttile(t);
            imagesc(ax, distance_km, squeeze(g.z), values, AlphaData=~isnan(values))
            hold(ax, "on")
            beam_range = linspace(0, 250e3, 200);
            beam_height = sqrt(beam_range .^ 2 + effective_radius ^ 2 + 2 * beam_range * effective_radius * sind(0.5)) - effective_radius + radar_selected.height(k);
            plot(ax, beam_range / 1000, beam_height, "w--", -beam_range / 1000, beam_height, "w--", LineWidth=0.8)
            set(ax, YDir="normal", Color=[0.85 0.85 0.85], FontSize=12, XLim=[-200 200], YLim=[0 3000])
            if i_row == 3
                clim(ax, [-1 1]); colormap(ax, signed_map)
                xlabel(ax, "Distance from radar, W (−) to E (+) (km)")
            else
                clim(ax, [0 150]); colormap(ax, parula)
            end
            if i_row == 1, title(ax, section_radars(i_k) + sprintf(" (%.2f°N)", g.lat(y_id)), FontSize=14), end
            if i_k == 1, ylabel(ax, "Height (m a.s.l.)"), end
            if i_k == numel(section_radars)
                cb = colorbar(ax);
                cb.FontSize = 11;
                if i_row == 3, cb.Label.String = "log_2 ratio"; else, cb.Label.String = "Density (birds km^{-3})"; end
            end
        end
    end
    drawnow
    frame_image = getframe(fig);
    [indexed, color_map] = rgb2ind(frame_image.cdata, 256);
    if i_t == 1
        imwrite(indexed, color_map, gif_file, LoopCount=Inf, DelayTime=1.2)
    else
        imwrite(indexed, color_map, gif_file, WriteMode="append", DelayTime=1.2)
    end
    if any(i_t == [3 6])
        bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "cfp_error_sections_" + string(c.times(i_t), "HHmm") + ".png"), Resolution=150)
    end
end
close(fig)

%% Local functions ----------------------------------------------------------
function prior_active = sandbox_prior(S)
% The radar-profile prior of the production inversion
inv = S.c.inversion;
evalc("prior_active = bmm4d.prior.build_density_prior_fast(S.g, S.radar_scans, S.observations, " + ...
    """lambda_profile"", inv.lambda_profile, ""lambda_time_profile"", inv.lambda_time_profile, " + ...
    """maxit"", inv.maxit, ""active_id"", S.active_id, ""interpolation"", inv.prior_interpolation, " + ...
    """length_scale_m"", inv.prior_length_scale_m, ""include_water"", true);");
end
