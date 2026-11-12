%% Sandbox: are the post-estimator filters of clean_volume still needed?
% clean_volume applied two filters to the density sweeps of read_rd before
% thinning (removed on 3 November 2026 after this check): a density cap
% (3000 birds km^-3) and a patch filter that removed connected patches of
% data > 0 smaller than 25 km^2. Both were designed for the
% paper's chain, where no-echo gates were 0. The gate estimator gives a
% positive density to every gate with information (grid from 0.1 birds
% km^-3), so data > 0 should equal ~isnan(data) and the patch filter can
% only remove valid islands between missing gates.
% Check on a few mid-night volumes of 12 May 2021, read serially at full
% resolution with the production options: whether data > 0 equals
% ~isnan(data), the minimum valid density, and the gates and share of
% birds (density x gate volume) removed by each filter.

%% Load packages and paths ------------------------------------------------
project_root = setup_paths();
addpath(fullfile(project_root, "examples"))
c = case_settings("night_20210512");
manifest = load(fullfile(c.data_dir, "manifest.mat"), "raw_files", "radars", "times");
check_radars = ["KDMX" "KLSX" "KMPX" "KICT"];
i_t = round(numel(manifest.times) / 2);

%% Clean each volume and apply the two filters by hand -------------------
check = table();
for radar = check_radars
    raw_file = manifest.raw_files(manifest.radars == radar, i_t);
    scan_cell = bmm4d.radar.read_rd(raw_file, c.read_options{:});
    sweeps = scan_cell{1}.dens.sweeps;
    positive_is_valid = true;
    min_valid = inf;
    max_valid = 0;
    n_cap = 0;
    n_patch = 0;
    bird_sum = 0;
    bird_cap = 0;
    bird_patch = 0;
    for i_elev = 1:numel(sweeps)
        data = sweeps(i_elev).data;
        range_m = sweeps(i_elev).range(:);
        az_step_rad = deg2rad(median(diff(sweeps(i_elev).az(:))));
        valid = ~isnan(data);
        positive_is_valid = positive_is_valid && isequal(data > 0, valid);
        min_valid = min([min_valid; data(valid)]);
        max_valid = max([max_valid; data(valid)]);
        % Birds per gate up to a constant (r^2 dAz dr): only shares are used
        bird = data .* range_m.^2 * az_step_rad * median(diff(range_m));
        bird_sum = bird_sum + sum(bird(valid));
        % Density cap, then patch filter, as in clean_volume (filter_sweep)
        cap = data > 3000;
        data(cap) = nan;
        gate_area_km2 = range_m * az_step_rad * median(diff(range_m)) / 1e6;
        patches = bwconncomp(data > 0, 8);
        patch_area = cellfun(@(ids) sum(gate_area_km2(mod(ids - 1, numel(range_m)) + 1)), patches.PixelIdxList);
        small = vertcat(patches.PixelIdxList{patch_area < 25});
        n_cap = n_cap + nnz(cap);
        n_patch = n_patch + numel(small);
        bird_cap = bird_cap + sum(bird(cap));
        bird_patch = bird_patch + sum(bird(small));
    end
    check = [check; table(radar, manifest.times(i_t), positive_is_valid, min_valid, max_valid, ...
        n_cap, bird_cap / bird_sum, n_patch, bird_patch / bird_sum, ...
        VariableNames=["radar" "time" "positive_is_valid" "min_density" "max_density" ...
        "n_cap" "cap_bird_share" "n_patch" "patch_bird_share"])]; %#ok<AGROW>
end
disp(check)
