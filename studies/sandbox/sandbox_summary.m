function summary = sandbox_summary(results, test_name)
%SANDBOX_SUMMARY Table and comparison figure for a set of sandbox runs.
%   RESULTS is a cell array of sandbox_run outputs. Writes
%   studies/sandbox/figures/<test_name>.csv and <test_name>.png: the
%   radar-distance curves (mean over frames), the bias per frame, and the
%   trade-off between radar bias and trusted misfit. Success is a bias near
%   0 (flat curve at 1) without a higher trusted misfit or a lower detail.

summary = table(cellfun(@(r) r.name, results), cellfun(@(r) r.bias_200km, results), ...
    cellfun(@(r) r.misfit, results), cellfun(@(r) r.misfit_trusted, results), cellfun(@(r) r.detail, results), ...
    cellfun(@(r) r.negative_share, results), cellfun(@(r) r.seconds, results), cellfun(@(r) r.iterations, results), ...
    VariableNames=["option" "bias_200km" "misfit" "misfit_trusted" "detail" "negative_share" "seconds" "iterations"]);
disp(summary)
figure_dir = fullfile(setup_paths(), "studies", "sandbox", "figures");
writetable(summary, fullfile(figure_dir, test_name + ".csv"))

fig = create_export_figure(3600, 1200, DefaultAxesFontName="Helvetica", DefaultTextFontName="Helvetica");
t = tiledlayout(fig, 1, 3, TileSpacing="tight", Padding="tight");
title(t, "Sandbox: " + strrep(test_name, "_", " "), FontSize=14, Color="k")
colors = lines(numel(results));
center = (results{1}.edges_km(1:end-1) + results{1}.edges_km(2:end)) / 2;
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_r = 1:numel(results)
    plot(ax, center, mean(results{i_r}.ratio, 2, "omitnan"), "-o", Color=colors(i_r, :), LineWidth=2, ...
        DisplayName=results{i_r}.name)
end
yline(ax, 1, "k--", HandleVisibility="off"); ylim(ax, [0.7 1.3])
xlabel(ax, "Distance to the nearest radar (km)"); ylabel(ax, "Column density / regional mean")
legend(ax, Location="southeast", FontSize=10, Interpreter="none"); set(ax, FontSize=12)
title(ax, "Radar-distance curve, mean over frames", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_r = 1:numel(results)
    plot(ax, results{i_r}.bias_frame, "-o", Color=colors(i_r, :), LineWidth=2)
end
xlabel(ax, "Frame"); ylabel(ax, "Radar-distance bias (within 200 km)"); set(ax, FontSize=12)
title(ax, "Bias per frame", FontWeight="normal")
ax = nexttile(t); hold(ax, "on"); box(ax, "on"); grid(ax, "on")
for i_r = 1:numel(results)
    scatter(ax, results{i_r}.misfit_trusted, results{i_r}.bias_200km, 120 - 8 * i_r, colors(i_r, :), "filled", ...
        MarkerEdgeColor="k")
end
ylim(ax, [0 max(ylim(ax))])
xlabel(ax, "Trusted misfit (30-120 km rows)"); ylabel(ax, "Radar-distance bias"); set(ax, FontSize=12)
title(ax, "Trade-off (best: bottom left; colours as left)", FontWeight="normal")
export_fixed_figure(fig, fullfile(figure_dir, test_name + ".png"), Resolution=150)
close(fig)
end
