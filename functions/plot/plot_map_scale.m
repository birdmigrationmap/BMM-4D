function plot_map_scale(ax, start_xy, length_plot, label, line_color, background_color, font_size)
% Draw a map scale with end caps and a centered distance label.

hold(ax, "on")
tick_height = 0.009 .* diff(ylim(ax));
x = start_xy(1) + [0 length_plot];
y = start_xy(2);
bar_x = [x nan x(1) x(1) nan x(2) x(2)];
bar_y = [y y nan y-tick_height y+tick_height nan y-tick_height y+tick_height];
plot(ax, bar_x, bar_y, Color=background_color, LineWidth=2.5)
plot(ax, bar_x, bar_y, Color=line_color, LineWidth=1)
text(ax, mean(x), y + 1.8 .* tick_height, label, FontName="Helvetica", FontSize=font_size, ...
    Color=line_color, BackgroundColor=background_color, Margin=0.5, ...
    HorizontalAlignment="center", VerticalAlignment="bottom")
end
