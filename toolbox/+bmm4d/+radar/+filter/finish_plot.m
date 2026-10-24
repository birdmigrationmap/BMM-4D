function finish_plot(fig, parms, suffix)
%FINISH_PLOT Export a filter's diagnostic figure, or show it.
%   With "export" in parms.plot_opts the figure is saved as
%   <export_dir>/<filename><SUFFIX> and closed; otherwise it is shown.

if ismember("export", parms.plot_opts)
    bmm4d.plot.export_fixed_figure(fig, fullfile(parms.export_dir, parms.filename + suffix), Resolution=300)
    close(fig)
else
    set(fig, "Visible", "on")
end
end
