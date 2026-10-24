function plot_cleaning_case(case_row, method_names, read_options, z_edges, output_file)
%PLOT_CLEANING_CASE Review figure comparing rain-removal methods on one volume.
%   Row 1: raw reflectivity, radial velocity, rhoHV and ZDR at 0.5 deg.
%   Row 2: bird density at 0.5 deg after each method, and where they differ.
%   Row 3: where methods differ at about 1.5 and 3.5 deg, the density
%   profiles; screening metrics are in the title. Cleaning runs at full gate
%   resolution so the masks are visible.

max_range = 150000;
file = case_row.file;

%% Read raw fields and clean with every method ----------------------------
radar = rsl2mat_wrap(file, struct("max_elev", 7));
read_options = [read_options(:)' {"data_res", [0 0]}];
scans = cell(1, numel(method_names));
for i_m = 1:numel(method_names)
    scan_cell = read_rd(file, "removeRain", method_names(i_m), read_options{:});
    scans{i_m} = scan_cell{1};
end
sweep_elev = [scans{1}.dens.sweeps.elev];

fig = create_export_figure(3600, 2700, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
% Keep the top 7% of the canvas for the two-line case title.
t = tiledlayout(fig, 3, 4, TileSpacing="compact", Padding="compact", OuterPosition=[0 0 1 0.93]);
title_lines = [ ...
    sprintf("%s  |  %s  |  %s %s UTC  |  %s", case_row.case_id, strrep(case_row.group, "_", " "), case_row.radar, ...
        string(radar.datetime, "yyyy-MM-dd HH:mm"), case_row.season)
    sprintf("Signal removed: rho %.1f%%, MistNet %.1f%%  |  VID (birds km^{-2}): none %.0f, rho %.0f, MistNet %.0f", ...
        100 * case_row.removed_rho, 100 * case_row.removed_mistnet, case_row.vid_none, case_row.vid_rho, case_row.vid_mistnet)];
annotation(fig, "textbox", [0 0.93 1 0.07], String=title_lines, EdgeColor="none", ...
    HorizontalAlignment="center", VerticalAlignment="middle", FontSize=17, FontWeight="bold", FontName="Helvetica")

%% Row 1: raw fields ------------------------------------------------------
raw_fields = ["dz" "vr" "rh" "dr"];
raw_titles = ["Z (dBZ)" "V_r (m s^{-1})" "\rho_{HV} (%)" "Z_{DR} (dB)"];
for i_f = 1:numel(raw_fields)
    ax = nexttile(t);
    sweeps = unique_elev_sweeps(radar, raw_fields(i_f));
    [~, i_sweep] = min(abs([sweeps.elev] - 0.5));
    [data, sweep_range, sweep_az] = sweep2mat(sweeps(i_sweep));
    [cmap, c_lim] = colormap_wct(raw_fields(i_f));
    if raw_fields(i_f) == "rh"
        data = 100 .* data;  % WCT rhoHV colormap is in percent
    elseif raw_fields(i_f) == "vr"
        c_lim = [-30 30];
    end
    plot_ppi(ax, data, sweep_az, sweep_range, max_range, cmap, c_lim, ...
        raw_titles(i_f) + sprintf(" %.1f°", sweeps(i_sweep).elev));
end

%% Row 2: density after each method and their difference ----------------
density_ref = scans{1}.dens.sweeps(1).data;
c_lim = [0 max(quantile(density_ref(:), 0.99), eps)];
for i_m = 1:numel(method_names)
    ax = nexttile(t);
    sweep = scans{i_m}.dens.sweeps(1);
    plot_ppi(ax, sweep.data, sweep.az, sweep.range, max_range, parula(256), c_lim, ...
        "Density, " + method_names(i_m));
end
plot_mask_comparison(nexttile(t), scans, 1, max_range)

%% Row 3: higher elevations, profiles, metrics ----------------------------
for target_elev = [1.5 3.5]
    [~, i_elev] = min(abs(sweep_elev - target_elev));
    plot_mask_comparison(nexttile(t), scans, i_elev, max_range)
end

ax = nexttile(t, [1 2]);
hold(ax, "on")
z_mid = (z_edges(1:end - 1) + z_edges(2:end)) / 2;
line_colors = [0.5 0.5 0.5; 0.85 0.4 0.1; 0.15 0.45 0.8];
for i_m = 1:numel(method_names)
    plot(ax, scan_density_profile(scans{i_m}, z_edges), z_mid, LineWidth=1.5, ...
        Color=line_colors(i_m, :), DisplayName=method_names(i_m))
end
xlabel(ax, "Mean density (birds km^{-3})")
ylabel(ax, "Height above radar (m)")
title(ax, "Density profiles")
legend(ax, Location="northeast", Box="off")
grid(ax, "on")
style_axes(ax)


export_fixed_figure(fig, output_file, Resolution=300)
close(fig)
end

function plot_mask_comparison(ax, scans, i_elev, max_range)
% Categories: kept by both, removed by rho only, by MistNet only, by both.
% Gates without valid density before rain removal stay blank.
density_ref = scans{1}.dens.sweeps(i_elev).data;
removed_rho = isnan(scans{2}.dens.sweeps(i_elev).data);
removed_mistnet = isnan(scans{3}.dens.sweeps(i_elev).data);
category = 1 + removed_rho + 2 .* removed_mistnet;
category(isnan(density_ref)) = nan;
sweep = scans{1}.dens.sweeps(i_elev);
cb = plot_ppi(ax, category, sweep.az, sweep.range, max_range, ...
    [0.82 0.82 0.82; 0.85 0.4 0.1; 0.15 0.45 0.8; 0.2 0.2 0.2], [0.5 4.5], sprintf("Removed %.1f°", sweep.elev));
set(cb, Ticks=1:4, TickLabels=["kept" "rho" "MN" "both"])
end

function cb = plot_ppi(ax, data, sweep_az, sweep_range, max_range, cmap, c_lim, ttl)
% Map one sweep around the radar in km; blank where there is no data.
[Z, x, y] = mat2cart(data, sweep_az, sweep_range, 1000, max_range, "nearest");
imagesc(ax, x ./ 1000, y ./ 1000, Z, AlphaData=~isnan(Z))
axis(ax, "xy", "equal")
axis(ax, [-1 1 -1 1] .* max_range ./ 1000)
set(ax, XTick=-150:150:150, YTick=-150:150:150, Box="on")
colormap(ax, cmap)
clim(ax, c_lim)
cb = colorbar(ax, Color="k", FontSize=11);
title(ax, ttl, FontWeight="normal")
style_axes(ax)
end

function style_axes(ax)
set(ax, FontName="Helvetica", FontSize=13, Color="w")
ax.Title.FontSize = 14;
end
