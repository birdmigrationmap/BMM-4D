%% Measure compact mosaic MAT-file size and selective load time

project_root = setup_paths();
file_patterns = [ ...
    fullfile(project_root, "06_mosaic_synthetic", "data", "case_*.mat")
    fullfile(project_root, "07_realistic_synthetic", "data", "case_*.mat")
    fullfile(project_root, "10_mosaic_reconstruction", "data", "case.mat")];

for i_pattern = 1:numel(file_patterns)
    files = dir(file_patterns(i_pattern));
    if contains(file_patterns(i_pattern), "06_mosaic_synthetic") || contains(file_patterns(i_pattern), "07_realistic_synthetic")
        files = files(~cellfun("isempty", regexp({files.name}, '^case_\d{2}\.mat$', "once")));
    end
    for i_file = 1:numel(files)
        file = fullfile(files(i_file).folder, files(i_file).name);
        file_info = whos("-file", file);
        variables = string({file_info.name});
        fprintf("%s | %.1f MB | %s\n", files(i_file).name, files(i_file).bytes / 1024^2, strjoin(variables, ", "))

        inversion_variables = cellstr(intersect(["observations" "grid" "radar_scans"], variables, "stable"));
        load_tic = tic;
        loaded = load(file, inversion_variables{:});
        load_seconds = toc(load_tic);
        loaded_memory = whos("loaded");
        fprintf("  inversion inputs: %.2f s, %.1f MB retained\n", load_seconds, loaded_memory.bytes / 1024^2)
        clear loaded

        plot_variables = cellstr(intersect(["truth_sequences" "grid" "radar_scans"], variables, "stable"));
        load_tic = tic;
        loaded = load(file, plot_variables{:});
        load_seconds = toc(load_tic);
        loaded_memory = whos("loaded");
        fprintf("  plotting inputs:  %.2f s, %.1f MB retained\n", load_seconds, loaded_memory.bytes / 1024^2)
        clear loaded
    end
end
