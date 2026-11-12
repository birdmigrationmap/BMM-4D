%% Gate density J: a reference set of rain-free migration scans with CFP (2021-2025)
% Selects one NEXRAD volume per night, three nights per radar (sparse,
% medium and dense migration), downloads it, checks it and keeps or
% rejects it. Every file of 2021 on records the clutter filter (CFP).
% 1. Nights from the Dark Ecology daily time series (darkeco-dataset, daily
%    file in local_data_dir/gate_density/darkeco/daily): spring (15 Mar-31 May)
%    and autumn (15 Aug-15 Nov), rain on less than 1% of the night's
%    scattering volumes, less than 10% of the night missing. Night-mean bird
%    reflectivity in three strata (tertiles above the 10th percentile).
% 2. Radars: beam blockage of the lowest cut within 50 km at most 5% (90th
%    percentile), water within 50 km at most 10%, antenna at most 50 m above
%    the ground (BMM-within data/reference).
% 3. Per radar and stratum, nights in random order until one passes: the
%    volume closest to 3 h after local sunset is downloaded and kept if the
%    lowest cut has CFP, has no rain (at most 0.2% of gates within 5-80 km
%    with REF >= 20 dBZ and RHO >= 0.97 at cuts below 1.6 deg) and has birds
%    (at least 5 birds km^-3 on average at 20-60 km). Rejected files are
%    deleted at once. Progress is saved after every night, so the script can
%    be stopped and run again.
% Outputs: local_data_dir/gate_density/post2020_scans (volumes, scans.mat,
% quicklook/*.png) and studies/gate_density/figures/GD_J_selected_scans.csv.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "studies", "gate_density", "functions"))
figure_dir = fullfile(project_root, "studies", "gate_density", "figures");

%% Settings ---------------------------------------------------------------
reference_dir = "/Users/rafnuss/Library/CloudStorage/OneDrive-Vogelwarte/BMM-within/data/reference";
daily_dir = fullfile(bmm4d.util.local_data_dir(), "gate_density", "darkeco", "daily");
out_dir = fullfile(bmm4d.util.local_data_dir(), "gate_density", "post2020_scans");
quicklook_dir = fullfile(out_dir, "quicklook");
years = 2021:2025;
n_stratum = 3;
max_attempts = 8;                          % nights tried per radar and stratum
hours_after_sunset = 3;
max_blockage = 0.05;
max_water = 0.10;
max_antenna_agl_m = 50;
max_rain_share = 0.002;
min_bird_density = 5;                      % birds km^-3 at 20-60 km, lowest cut
max_range_m = 80000;
if ~isfolder(quicklook_dir)
    mkdir(quicklook_dir)
end

%% Radars -------------------------------------------------------------------
radar = readtable(fullfile(reference_dir, "radar.csv"), TextType="string");
blockage = load(fullfile(reference_dir, "beam_blockage", "beam_blockage"));
blockage_q90 = quantile(reshape(blockage.bb(:, blockage.bb_range < 50000, :), [], size(blockage.bb, 3)), 0.9)';
radar.antenna_agl = radar.height - radar.dem;
radar = radar(blockage_q90 <= max_blockage & radar.onwater50 <= max_water & radar.antenna_agl <= max_antenna_agl_m, :);

%% Candidate nights ---------------------------------------------------------
nights = table();
for y = years
    daily = readtable(fullfile(daily_dir, y + "-daily.csv"), TextType="string");
    nights = [nights; daily(daily.period == "night", :)]; %#ok<AGROW>
end
day_of_year = day(nights.date, "dayofyear");
nights.season = repmat("", height(nights), 1);
nights.season(day_of_year >= 74 & day_of_year <= 151) = "spring";
nights.season(day_of_year >= 227 & day_of_year <= 319) = "autumn";
nights.mean_refl = nights.reflectivity_hours ./ nights.period_length;     % cm^2 km^-2
nights = nights(nights.season ~= "" & nights.fraction_rain < 0.01 & nights.fraction_missing < 0.1 & ...
    ismember(nights.station, radar.name) & nights.mean_refl > 0, :);
nights = nights(nights.mean_refl >= prctile(nights.mean_refl, 10), :);
nights.stratum = discretize(nights.mean_refl, [0 prctile(nights.mean_refl, [100 / 3 200 / 3]) Inf]);
fprintf("[select] %d radars, %d candidate nights\n", height(radar), height(nights))

%% Download, check, keep or reject -----------------------------------------------
% Resume from earlier runs
progress_file = fullfile(out_dir, "scans.mat");
if isfile(progress_file)
    load(progress_file, "scans", "tried")
else
    scans = table();
    tried = strings(0, 1);
end
for i_r = 1:height(radar)
    station = radar.name(i_r);
    for k = 1:n_stratum
        if ~isempty(scans) && any(scans.station == station & scans.stratum == k)
            continue
        end
        candidate = nights(nights.station == station & nights.stratum == k, :);
        candidate = candidate(~ismember(station + string(candidate.date, "yyyyMMdd"), tried), :);
        rng(100 * i_r + k)
        candidate = candidate(randperm(height(candidate), min(height(candidate), max_attempts)), :);
        for i_c = 1:height(candidate)
            night = candidate(i_c, :);
            tried(end + 1, 1) = station + string(night.date, "yyyyMMdd"); %#ok<SAGROW>
            % Night starting on the local date: the sunset and sunrise around local midnight, in UTC
            [sunset, sunrise] = bmm4d.util.getSunriseSunset(radar.lat(i_r), radar.lon(i_r), ...
                datetime(night.date) + hours(24 - radar.lon(i_r) / 15));
            target_time = sunset + hours(hours_after_sunset);
            if target_time > sunrise - hours(1)
                continue
            end
            try
                raw_file = gd_download_level2(station, target_time, out_dir);
            catch err
                fprintf("[select] %s %s download failed: %s\n", station, string(night.date), err.message)
                continue
            end
            if raw_file == ""
                continue
            end
            % Checks on the low cuts
            verdict = "keep";
            try
                level2 = bmm4d.radar.read_level2_moments(raw_file, ["REF" "RHO" "CFP"], 1.6);
                level2 = level2(arrayfun(@(s) ~isempty(s.REF), level2));
                if isempty(level2(1).CFP)
                    verdict = "no CFP";
                else
                    n_rain = 0;
                    n_gate = 0;
                    for i_e = 1:numel(level2)
                        if isempty(level2(i_e).RHO) || size(level2(i_e).RHO, 2) ~= size(level2(i_e).REF, 2)
                            continue
                        end
                        n_common = min(size(level2(i_e).REF, 1), size(level2(i_e).RHO, 1));
                        in_range = level2(i_e).REF_range(1:n_common) >= 5000 & level2(i_e).REF_range(1:n_common) <= max_range_m;
                        rain = level2(i_e).REF(in_range, :) >= 20 & level2(i_e).RHO(in_range, :) >= 0.97;
                        n_rain = n_rain + nnz(rain);
                        n_gate = n_gate + numel(rain);
                    end
                    density = z_to_refl(idb(level2(1).REF)) / 11;
                    density(isnan(density)) = 0;
                    density = min(density, 3000);
                    far = level2(1).REF_range >= 20000 & level2(1).REF_range <= 60000;
                    bird_density = mean(density(far, :), "all");
                    if n_gate == 0 || n_rain / n_gate > max_rain_share
                        verdict = "rain";
                    elseif bird_density < min_bird_density
                        verdict = "few birds";
                    end
                end
            catch err
                verdict = "unreadable: " + err.message;
            end
            fprintf("[select] %s %s stratum %d: %s\n", station, string(night.date), k, verdict)
            if verdict ~= "keep"
                delete(raw_file)
                continue
            end
            % Keep: record and quick-look of the lowest cut (reflectivity and RHO)
            info = dir(raw_file);
            time = datetime(extractBetween(info.name, station, "_V"), InputFormat="yyyyMMdd_HHmmss");
            scans = [scans; table(station, time, string(raw_file), radar.antenna_agl(i_r), night.date, night.season, ...
                k, night.mean_refl, bird_density, n_rain / n_gate, VariableNames=["station" "time" "raw_file" ...
                "antenna_agl" "night" "season" "stratum" "mean_refl" "bird_density" "rain_share"])]; %#ok<AGROW>
            fig = bmm4d.plot.create_export_figure(1200, 560, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
            t = tiledlayout(fig, 1, 2, TileSpacing="tight", Padding="tight");
            in_map = level2(1).REF_range <= max_range_m;
            [az_grid, r_grid] = meshgrid(0:360 / size(level2(1).REF, 2):360, level2(1).REF_range(in_map) / 1000);
            panels = {level2(1).REF(in_map, :), level2(1).RHO(1:min(nnz(in_map), size(level2(1).RHO, 1)), :)};
            panel_limits = [-10 30; 0.3 1];
            panel_names = ["REF (dBZ)" "RHO"];
            for i_p = 1:2
                ax = nexttile(t);
                values = panels{i_p};
                rows = 1:size(values, 1);
                pcolor(ax, r_grid(rows, :) .* sind(az_grid(rows, :)), r_grid(rows, :) .* cosd(az_grid(rows, :)), values(:, [1:end 1]))
                shading(ax, "flat"); axis(ax, "equal", "off"); colorbar(ax);
                colormap(ax, parula); clim(ax, panel_limits(i_p, :))
                title(ax, panel_names(i_p) + " " + station + " " + string(time, "yyyy-MM-dd HH:mm"), FontWeight="normal")
            end
            bmm4d.plot.export_fixed_figure(fig, fullfile(quicklook_dir, station + "_" + string(time, "yyyyMMdd_HHmm") + ".png"), Resolution=100)
            close(fig)
            break
        end
        save(progress_file, "scans", "tried")
    end
    fprintf("[select] %d/%d radars, %d scans kept\n", i_r, height(radar), height(scans))
end
writetable(scans, fullfile(figure_dir, "GD_J_selected_scans.csv"))
fprintf("[select] done: %d scans\n", height(scans))
