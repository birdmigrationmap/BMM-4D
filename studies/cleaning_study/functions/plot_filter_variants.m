function plot_filter_variants(raw_file, variant_files, variant_names, output_file)
%PLOT_FILTER_VARIANTS Compare cleaned volumes of one radar volume.
%   Rows: the lowest sweep and the sweep closest to 1.5 deg. Columns: raw
%   reflectivity, raw rhoHV, then the density of each cleaning variant
%   (grey: zero density from no-echo gates; white: missing). A last panel
%   shows the mean density profile of each variant.

raw = bmm4d.radar.rsl2mat_wrap(char(raw_file), struct("max_elev", 7));
variants = cellfun(@(f) load(f, "volume").volume, cellstr(variant_files));
n_var = numel(variants);
elev = [variants(1).dens.sweeps.elev];
[~, i_upper] = min(abs(elev - 1.5));
sweep_ids = unique([1 i_upper]);

fig = bmm4d.plot.create_export_figure(450 * (n_var + 3), 420 * numel(sweep_ids), Visible="off");
t = tiledlayout(fig, numel(sweep_ids), n_var + 3, TileSpacing="tight", Padding="tight");
density_cmap = [0.72 0.72 0.72; parula(255)];
for k = 1:numel(sweep_ids)
    target = elev(sweep_ids(k));
    raw_fields = ["dz" "rh"];
    raw_lims = {[-10 35], [0.3 1]};
    for f = 1:2
        sweeps = unique_elev_sweeps(raw, char(raw_fields(f)));
        [~, j] = min(abs([sweeps.elev] - target));
        [data, range_m, az] = sweep2mat(sweeps(j));
        ax = nexttile(t);
        polar_surface(ax, data, range_m, az);
        clim(ax, raw_lims{f}); colorbar(ax);
        title(ax, sprintf("raw %s %.1f°", raw_fields(f), sweeps(j).elev))
    end
    for i_v = 1:n_var
        sweep = variants(i_v).dens.sweeps(sweep_ids(k));
        data = fill_thinned(sweep.data);
        data(data == 0) = -1;   % grey
        ax = nexttile(t);
        polar_surface(ax, data, sweep.range, sweep.az);
        colormap(ax, density_cmap); clim(ax, [-1 150]); colorbar(ax);
        title(ax, variant_names(i_v))
    end
    ax = nexttile(t);
    if k == 1
        hold(ax, "on"); box(ax, "on")
        for i_v = 1:n_var
            [z, profile] = mean_profile(variants(i_v));
            plot(ax, profile, z / 1000, LineWidth=1.5, DisplayName=sprintf("%s (VID %.0f)", variant_names(i_v), sum(profile, "omitnan") * 0.2));
        end
        xlabel(ax, "mean density (birds km^{-3})"); ylabel(ax, "height above radar (km)")
        legend(ax, Location="northeast")
    else
        axis(ax, "off")
    end
end
title(t, sprintf("%s %s", variants(1).station, string(variants(1).datetime, "yyyy-MM-dd HH:mm")) + " UTC")
bmm4d.plot.export_fixed_figure(fig, output_file, Resolution=300)
close(fig)
end

function data = fill_thinned(data)
% Repeat each kept gate over its thinned neighbours so the map is readable
% (display only: a missing kept gate can be covered by its neighbour).
kept_rows = find(any(~isnan(data), 2));
for i = kept_rows'
    kept_cols = find(~isnan(data(i, :)));
    if numel(kept_cols) > 1
        data(i, :) = fillmissing(data(i, :), "previous", MaxGap=median(diff(kept_cols)));
    end
end
if numel(kept_rows) > 1
    data = fillmissing(data, "previous", 1, MaxGap=median(diff(kept_rows)));
end
end

function [z_center, profile] = mean_profile(volume)
% Mean density in 200 m bins of beam height, all sweeps together.
edges = 0:200:3600;
z_center = edges(1:end-1) + 100;
values = [];
heights = [];
for sweep = volume.dens.sweeps
    [~, h] = slant2ground(sweep.range(:), sweep.elev);
    h = repmat(h, 1, size(sweep.data, 2));
    ok = ~isnan(sweep.data);
    values = [values; sweep.data(ok)]; %#ok<AGROW>
    heights = [heights; h(ok)]; %#ok<AGROW>
end
bin = discretize(heights, edges);
profile = accumarray(bin(~isnan(bin)), values(~isnan(bin)), [numel(z_center) 1], @mean, nan)';
end

function polar_surface(ax, data, range_m, az)
% Map a polar sweep to x-y km within 150 km.
range_m = range_m(:);
keep = range_m <= 150000;
[A, R] = meshgrid(az(:)', range_m(keep));
surface(ax, R .* sind(A) / 1000, R .* cosd(A) / 1000, zeros(size(A)), data(keep, :), EdgeColor="none");
axis(ax, "equal"); xlim(ax, [-150 150]); ylim(ax, [-150 150]); box(ax, "on")
set(ax, Color="w", XTick=[-100 0 100], YTick=[-100 0 100])
end
