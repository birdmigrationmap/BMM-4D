function gd_plot_band(ax, x, values, color, name)
%GD_PLOT_BAND Mean line and min-max band of replicate curves (columns).
fill(ax, [x; flipud(x)], [min(values, [], 2); flipud(max(values, [], 2))], color, FaceAlpha=0.15, EdgeColor="none", ...
    HandleVisibility="off")
plot(ax, x, mean(values, 2), "-o", Color=color, LineWidth=2, DisplayName=name)
end
