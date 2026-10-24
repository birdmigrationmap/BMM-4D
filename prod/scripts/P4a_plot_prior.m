%% Production test 4a: density prior for review
% Builds the density prior of the case with the paper interpolation
% (nearest-radar IDW) and with the smooth Gaussian weights, and plots the
% column density at six times and a west-east vertical section, so the
% prior can be checked before the inversion. Needs P3 (case.mat).

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "prod", "functions"))
if ~exist("case_name", "var")
    case_name = "keywest";
end
c = prod_case(case_name);
inv = c.inversion;
load(fullfile(c.data_dir, "case.mat"), "observations", "g", "radar_scans", "radar_selected")
n_time = size(radar_scans, 2);
column_mask = [];
if ~isempty(c.domain)
    column_mask = build_domain_mask(g, c.domain);
end
active_mask_3d = build_active_mask(g, radar_scans, c.active_range_m, c.flight_layer_agl_m, column_mask);
active_id = find(active_mask_3d(:));

%% Prior with both interpolations ----------------------------------------
common = {"lambda_profile", inv.lambda_profile, "lambda_time_profile", inv.lambda_time_profile, ...
    "idw_neighbors", inv.density_prior_idw_neighbors, "idw_power", inv.density_prior_idw_power, ...
    "maxit", inv.maxit, "active_id", active_id};
variants = ["idw" "gaussian"];
variant_titles = ["Paper prior (IDW, 3 nearest radars, land fit)", ...
    sprintf("Smooth prior (Gaussian, L = %.0f km, land and water fit)", inv.prior_length_scale_m / 1000)];
prior_column = cell(1, 2);
prior_section = cell(1, 2);
[~, section_y] = min(abs(g.lat - mean(radar_selected.lat)));   % west-east section at the radars' mean latitude
for k = 1:2
    stage_tic = tic;
    prior_active = build_density_prior_fast(g, radar_scans, observations, common{:}, ...
        interpolation=variants(k), length_scale_m=inv.prior_length_scale_m, include_water=variants(k) == "gaussian");
    stage_log(c, "P4a_prior_" + variants(k), stage_tic)
    column = nan([g.sz(1:2) n_time], "single");
    section = nan([g.sz(1) g.sz(3) n_time], "single");
    for i_t = 1:n_time
        field = nan(g.sz, "single");
        field(active_id) = prior_active(:, i_t);
        column(:, :, i_t) = sum(field, 3, "omitnan") .* g.dz ./ 1000;
        section(:, :, i_t) = squeeze(field(:, section_y, :));
    end
    column(repmat(~any(active_mask_3d, 3), 1, 1, n_time)) = nan;
    prior_column{k} = column;
    prior_section{k} = section;
end
clear prior_active field

%% Plot column density and section ---------------------------------------
if ~isfolder(c.figure_dir)
    mkdir(c.figure_dir)
end
frame_ids = unique(round(linspace(1, n_time, 5)));
frame_ids = frame_ids(2:4);   % early, middle and late night
c_lim = [0 quantile(prior_column{2}(:), 0.995)];
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
land_color = [0.96 0.96 0.94];
water_color = [0.85 0.92 0.97];
background = permute(reshape(land_color, 1, 1, 3) .* ~g.mask_water + reshape(water_color, 1, 1, 3) .* g.mask_water, [2 1 3]);
mid = frame_ids(2);
fig = create_export_figure(3600, 3600, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 4, 3, TileSpacing="tight", Padding="tight");
for k = 1:2
    for i_f = 1:numel(frame_ids)
        ax = nexttile(t);
        plot_night_map(ax, g, background, prior_column{k}(:, :, frame_ids(i_f)), c_lim, states, radar_selected, c.figure_bounds, ones(g.sz(1:2)));
        yline(ax, g.lat(section_y), "w--", LineWidth=1)
        title(ax, variant_titles(k) + newline + string(c.times(frame_ids(i_f)), "dd MMM HH:mm") + " UTC", FontSize=12)
    end
end
cb = colorbar(ax, Color="k");
cb.Layout.Tile = "east";
cb.Label.String = "Prior column density (birds km^{-2})";
% Vertical sections (height above sea level vs longitude) at the middle frame
section_lim = [0 quantile(prior_section{2}(:, :, mid), 0.99, "all")];
for k = 1:2
    ax = nexttile(t, [1 3]);
    imagesc(ax, g.lon, squeeze(g.z), prior_section{k}(:, :, mid)', AlphaData=~isnan(prior_section{k}(:, :, mid)'))
    set(ax, YDir="normal", Color=[0.9 0.9 0.9], FontSize=12)
    colormap(ax, parula(256)); clim(ax, section_lim)
    scb = colorbar(ax, Color="k"); scb.Label.String = "Prior density (birds km^{-3})";
    xlim(ax, c.figure_bounds(1:2))
    title(ax, variant_titles(k) + sprintf(", west-east section at %.1f°N, ", g.lat(section_y)) + string(c.times(mid), "HH:mm") + " UTC", FontSize=12)
    xlabel(ax, "Longitude"); ylabel(ax, "Height above sea level (m)")
end
export_fixed_figure(fig, fullfile(c.figure_dir, "prior_review.png"), Resolution=300)
close(fig)
fprintf("[prod/%s] prior review in %s\n", c.name, fullfile(c.figure_dir, "prior_review.png"))
