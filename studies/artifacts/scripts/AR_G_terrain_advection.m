%% Artefacts G: advection over terrain (dev, follow-up of fix 2)
% Does the advection carry birds over the Appalachians without losing them
% into the ground? Three vertical rules for the advection links (ar_steps),
% each with its matching prior (ar_profile_prior):
% - ASL (production): a link keeps its height above sea level; a link
%   whose destination is inside a ridge is dropped, and the lee gets no
%   inflow from below the ridge top;
% - AGL: a link keeps its height above the ground; the column rises over a
%   ridge, but the layers pushed above the 3000 m grid top are dropped;
% - sigma: a link keeps its relative height between the ground and the top
%   of its column and carries the ratio of the columns' air thickness
%   (ar_solve_advection_lean), so a column squeezed over a ridge keeps its
%   birds: nothing goes into the ground or out of the top.
% Case ar_appalachians_bb (16 radars, 04:30-07:30 UTC, beam-blockage fix),
% production inversion settings with the white prior.
% Scores: links kept over terrain, ridge / valley column, column ratio per
% km of ground, misfit, and the prediction error of held-out mountain
% radars (each radar's rows get weight 0 and are predicted from the
% others: the misfit cannot rank the rules, a held-out radar can).
% Outputs: figures/AR_G_links.csv, AR_G_terrain_advection.csv,
% AR_G_holdout.csv, AR_G_terrain_advection.png, AR_G_section.png.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"), fullfile(project_root, "studies", "artifacts", "functions"))
figure_dir = fullfile(project_root, "studies", "artifacts", "figures");
times = datetime("2026-10-22 04:30"):minutes(30):datetime("2026-10-22 07:30");
c = ar_test_case("ar_appalachians_bb", [-89 -78 33.5 40], times);
S = ar_load(c);
inv = c.inversion;
dt_seconds = seconds(diff(c.times));
modes = ["asl" "agl" "sigma"];

%% Links and priors of each vertical rule ---------------------------------------
steps = struct();
priors = struct();
for mode = modes
    steps.(mode) = ar_steps(S.g, S.radar_scans, dt_seconds, S.active_mask_3d, inv, mode);
    priors.(mode) = ar_profile_prior(S.g, S.radar_scans, S.radar_profiles, S.active_mask_3d, inv.prior_length_scale_m, mode);
end

% Terrain crossings, from the sigma links (every source cell has one):
% uphill sources (destination ground > 300 m higher) and lee cells (in a
% column reached from ground > 300 m higher). Per rule: the share of
% uphill sources that keep a link and of lee cells that receive one.
dem = max(S.g.dem, 0);
nxy = prod(S.g.sz(1:2));
column_of = @(id) mod(double(id) - 1, nxy) + 1;
active_cell = false(S.g.sz);
active_cell(S.active_id) = true;
counts = zeros(numel(modes), 4);
for i_t = 1:numel(steps.sigma)
    s = steps.sigma{i_t};
    rise = dem(column_of(s.dest_id_valid)) - dem(column_of(s.source_id_valid));
    uphill = unique(s.source_id_valid(rise > 300));
    lee_column = false(S.g.sz(1:2));
    lee_column(unique(column_of(s.dest_id_valid(rise < -300)))) = true;
    lee = find(active_cell & repmat(lee_column, 1, 1, S.g.sz(3)));
    for i_m = 1:numel(modes)
        m = steps.(modes(i_m)){i_t};
        counts(i_m, :) = counts(i_m, :) + [nnz(ismember(uphill, m.source_id_valid)) numel(uphill) ...
            nnz(ismember(lee, m.dest_id_valid)) numel(lee)];
    end
end
links = table(modes', counts(:, 1) ./ counts(:, 2), counts(:, 3) ./ counts(:, 4), ...
    VariableNames=["rule" "uphill_sources_linked" "lee_cells_reached"]);
disp(links)
writetable(links, fullfile(figure_dir, "AR_G_links.csv"))

%% Solves ------------------------------------------------------------------------
variants = struct("name", {"ASL (production)", "AGL", "sigma", "sigma links, AGL prior"}, ...
    "steps", {"asl", "agl", "sigma", "sigma"}, "prior", {"asl", "agl", "sigma", "agl"});
results = cell(numel(variants), 1);
scores = table();
for i_v = 1:numel(variants)
    v = variants(i_v);
    results{i_v} = ar_solve(S, v.name, steps=steps.(v.steps), prior=priors.(v.prior));
    score = ar_terrain_score(S, results{i_v}.column);
    scores = [scores; table(string(v.name), score.ridge_valley, score.per_km, results{i_v}.misfit, ...
        results{i_v}.misfit_trusted, results{i_v}.seconds, results{i_v}.iterations, ...
        VariableNames=["solve" "ridge_valley" "per_km" "misfit" "misfit_30_120km" "seconds" "iterations"])]; %#ok<AGROW>
end
disp(scores)
writetable(scores, fullfile(figure_dir, "AR_G_terrain_advection.csv"))

%% Held-out mountain radars ----------------------------------------------------------
% Each radar's rows get weight 0 and the prior is rebuilt from the other
% radars' profiles; its rows at 30-120 km are then predicted from the other
% radars: relative L2 error and median |log ratio|.
held_out = intersect(["KMRX" "KGSP" "KFCX" "KRLX" "KJKL" "KHTX"], string(S.radar_selected.name), "stable");
holdout = table();
for radar_name = held_out(:)'
    for i_v = 1:numel(variants)
        v = variants(i_v);
        score = ar_holdout_score(S, radar_name, @(S_out, keep) ar_solve(S_out, v.name + " without " + radar_name, ...
            steps=steps.(v.steps), prior=ar_profile_prior(S.g, S.radar_scans(keep, :), S.radar_profiles(:, keep, :), ...
            S.active_mask_3d, inv.prior_length_scale_m, v.prior)));
        holdout = [holdout; table(radar_name, string(v.name), score.relative_error, score.median_abs_log_ratio, ...
            VariableNames=["radar" "solve" "relative_error" "median_abs_log_ratio"])]; %#ok<AGROW>
    end
end
disp(holdout)
disp(groupsummary(holdout, "solve", "mean", ["relative_error" "median_abs_log_ratio"]))
writetable(holdout, fullfile(figure_dir, "AR_G_holdout.csv"))

%% Maps ----------------------------------------------------------------------------
ar_plot_maps(S, results, [1 4 7], fullfile(figure_dir, "AR_G_terrain_advection.png"), ...
    bounds=[-88 -79 34 38.5], extra=dem, extra_title="Ground (m)", extra_limits=[0 1600]);

%% Section across the Smokies at 06:00 UTC -----------------------------------------
% West-east section at 35.6 N (Clingmans Dome 35.56 N, -83.50 E): density
% by height above sea level; the ground in gray. imagesc wants [z x], so
% the [x z] slice is transposed here.
[~, y_id] = min(abs(S.g.lat - 35.6));
x_id = find(S.g.lon >= -86 & S.g.lon <= -81);
i_t = find(c.times == datetime("2026-10-22 06:00"));
fig = bmm4d.plot.create_export_figure(2400, 300 + 450 * numel(results), DefaultAxesFontName="Helvetica", ...
    DefaultTextFontName="Helvetica");
layout = tiledlayout(fig, numel(results), 1, TileSpacing="tight", Padding="tight");
sections = cell(numel(results), 1);
for i_v = 1:numel(results)
    frame = nan(S.g.sz);
    frame(S.active_id) = results{i_v}.field(:, i_t);
    sections{i_v} = squeeze(frame(x_id, y_id, :));
end
limit = quantile(cell2mat(cellfun(@(s) s(:), sections, UniformOutput=false)), 0.99);
for i_v = 1:numel(results)
    ax = nexttile(layout);
    imagesc(ax, S.g.lon(x_id), S.g.z(:) / 1000, sections{i_v}', AlphaData=~isnan(sections{i_v}'));
    hold(ax, "on")
    ground = dem(x_id, y_id) / 1000;
    area(ax, S.g.lon(x_id), ground, FaceColor=[0.55 0.55 0.55], EdgeColor="none")
    set(ax, YDir="normal", Color="w", FontSize=13)
    colormap(ax, parula)
    clim(ax, [0 limit])
    ylim(ax, [0 3])
    ylabel(ax, "Height (km asl)")
    title(ax, string(results{i_v}.name) + ", 06:00 UTC, 35.6° N", FontSize=14, FontWeight="bold")
    cb = colorbar(ax);
    cb.Label.String = "Birds km^{-3}";
end
xlabel(ax, "Longitude (°)")
bmm4d.plot.export_fixed_figure(fig, fullfile(figure_dir, "AR_G_section.png"), Resolution=300);
close(fig)
