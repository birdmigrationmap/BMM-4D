%% Sandbox: production reconstruction before and after the observation weights
% Both from the production pipeline on the same observations (the sandbox
% night, volumes cleaned with the estimator sd): equal weights (before,
% prod/sandbox_equal_weights) and weights from the estimator's sd at the
% prior's prediction with a 1.5 birds km^-3 floor (after, prod/sandbox,
% bmm4d.solve.observation_weights).
% Rows: before, after, log2(after / before). Columns: column density from
% above (white: section lines), west-east sections through KDMX and KLSX,
% south-north through KDVN. One frame per time (gif), the middle frame
% also as PNG.
% Outputs: figures/production_weights.gif, production_weights.png.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"))
figure_dir = fullfile(project_root, "studies", "sandbox", "figures");
c = case_settings("sandbox");
load(fullfile(c.data_dir, "case.mat"), "g", "radar_selected")
after = load(fullfile(c.data_dir, "inversion.mat"), "active_field", "active_id");
before = load(fullfile(fileparts(c.data_dir), "sandbox_equal_weights", "inversion.mat"), "active_field", "active_id");
assert(isequal(after.active_id, before.active_id))
active_id = after.active_id;
fields = {before.active_field, after.active_field};
row_names = ["before (equal weights)" "after (sd weights)" "log_2 after / before"];
n_time = size(fields{1}, 2);

%% Geometry of the panels -------------------------------------------------------
[ix, iy, ~] = ind2sub(g.sz, active_id);
column_id = sub2ind(g.sz(1:2), ix, iy);
air = bmm4d.grid.air_fraction(g, active_id);
to_column = @(x) accumarray(column_id, double(x) .* air, [prod(g.sz(1:2)) 1]) * g.dz / 1000;   % birds km^-2
active_column = to_column(ones(numel(active_id), 1)) > 0;
states = shaperead(which("usastatehi.shp"), UseGeoCoords=true);
background = permute(reshape([0.96 0.96 0.94], 1, 1, 3) .* ~g.mask_water + reshape([0.85 0.92 0.97], 1, 1, 3) .* g.mask_water, [2 1 3]);
signed_map = interp1([0 0.5 1], [0.13 0.4 0.75; 1 1 1; 0.8 0.2 0.15], linspace(0, 1, 256));
height_m = squeeze(g.z);
sections = table(["KDMX" "KLSX" "KDVN"]', ["we" "we" "sn"]', VariableNames=["radar" "direction"]);
for i_k = 1:height(sections)
    k = find(radar_selected.name == sections.radar(i_k), 1);
    sections.lon(i_k) = radar_selected.lon(k);
    sections.lat(i_k) = radar_selected.lat(k);
    [~, sections.x_id(i_k)] = min(abs(g.lon - sections.lon(i_k)));
    [~, sections.y_id(i_k)] = min(abs(g.lat - sections.lat(i_k)));
end

%% Animation ------------------------------------------------------------------
file_stem = fullfile(figure_dir, "production_weights");
fig = bmm4d.plot.create_export_figure(3600, 2500, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
for i_t = 1:n_time
    clf(fig)
    t = tiledlayout(fig, 3, 1 + height(sections), TileSpacing="tight", Padding="tight");
    title(t, "Sandbox night, production run, " + string(c.times(i_t), "dd MMM yyyy HH:mm") + " UTC", FontSize=18, FontWeight="bold")
    volumes = cell(1, 2);
    for i_f = 1:2
        volumes{i_f} = nan(g.sz, "single");
        volumes{i_f}(active_id) = fields{i_f}(:, i_t);
    end
    for i_row = 1:3
        is_ratio = i_row == 3;
        % From above
        ax = nexttile(t);
        column = nan(g.sz(1:2));
        if is_ratio
            ratio = log2(max(to_column(fields{2}(:, i_t)), 1) ./ max(to_column(fields{1}(:, i_t)), 1));
            column(active_column) = ratio(active_column);
            bmm4d.plot.plot_night_map(ax, g, background, column, [-0.5 0.5], states, radar_selected, c.figure_bounds, ones(g.sz(1:2)));
            colormap(ax, signed_map)
        else
            values = to_column(fields{i_row}(:, i_t));
            column(active_column) = log10(max(values(active_column), 1));
            bmm4d.plot.plot_night_map(ax, g, background, column, log10([5 300]), states, radar_selected, c.figure_bounds, ones(g.sz(1:2)));
            colormap(ax, parula)
        end
        hold(ax, "on")
        for i_k = 1:height(sections)
            if sections.direction(i_k) == "we"
                plot(ax, sections.lon(i_k) + [-3 3], g.lat(sections.y_id(i_k)) * [1 1], "k-", LineWidth=1.5)
            else
                plot(ax, g.lon(sections.x_id(i_k)) * [1 1], sections.lat(i_k) + [-3 3], "k-", LineWidth=1.5)
            end
        end
        ylabel(ax, row_names(i_row), FontWeight="bold", FontSize=14)
        if i_row == 1, title(ax, "Column, 5-300 km^{-2}", FontSize=13), end
        if i_row == 3, title(ax, "Column, \pm0.5", FontSize=13), end
        % Sections
        for i_k = 1:height(sections)
            ax = nexttile(t);
            if is_ratio
                field = log2(max(volumes{2}, 1) ./ max(volumes{1}, 1));
            else
                field = volumes{i_row};
            end
            if sections.direction(i_k) == "we"
                values = squeeze(field(:, sections.y_id(i_k), :))';
                imagesc(ax, g.lon, height_m, values, AlphaData=~isnan(values))
                xlim(ax, sections.lon(i_k) + [-3 3])
                header = sprintf("W-E through %s (%.2f°N)", sections.radar(i_k), g.lat(sections.y_id(i_k)));
                x_label = "Longitude (°)";
            else
                values = squeeze(field(sections.x_id(i_k), :, :))';
                imagesc(ax, g.lat, height_m, values, AlphaData=~isnan(values))
                xlim(ax, sections.lat(i_k) + [-3 3])
                header = sprintf("S-N through %s (%.2f°W)", sections.radar(i_k), -g.lon(sections.x_id(i_k)));
                x_label = "Latitude (°)";
            end
            hold(ax, "on")
            xline(ax, ternary_position(sections, i_k), "k:", LineWidth=1.5)
            set(ax, YDir="normal", Color=[0.85 0.85 0.85], FontSize=12); ylim(ax, [0 3000])
            if is_ratio
                clim(ax, [-1 1]); colormap(ax, signed_map)
            else
                clim(ax, [0 150]); colormap(ax, parula)
            end
            if i_row == 1, title(ax, header, FontSize=14), end
            if i_k == 1, ylabel(ax, "Height (m a.s.l.)"), end
            if i_row == 3, xlabel(ax, x_label), end
            if i_k == height(sections)
                cb = colorbar(ax); cb.FontSize = 11;
                if is_ratio
                    cb.Label.String = "log_2 ratio";
                else
                    cb.Label.String = "Density (birds km^{-3})";
                end
            end
        end
    end
    drawnow
    frame_image = getframe(fig);
    [indexed, color_map] = rgb2ind(frame_image.cdata, 256);
    if i_t == 1
        imwrite(indexed, color_map, file_stem + ".gif", LoopCount=Inf, DelayTime=1.2)
    else
        imwrite(indexed, color_map, file_stem + ".gif", WriteMode="append", DelayTime=1.2)
    end
    if i_t == ceil(n_time / 2)
        bmm4d.plot.export_fixed_figure(fig, file_stem + ".png", Resolution=150)
    end
end
close(fig)

%% Local functions ----------------------------------------------------------
function position = ternary_position(sections, i_k)
% Radar position along a section
if sections.direction(i_k) == "we"
    position = sections.lon(i_k);
else
    position = sections.lat(i_k);
end
end
