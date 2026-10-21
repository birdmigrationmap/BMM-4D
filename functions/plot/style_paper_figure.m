function style_paper_figure(fig, font_size)
% Apply the common paper typography without changing scientific colors or limits.
set(findall(fig, "-property", "FontName"), FontName="Helvetica")
set(findall(fig, "-property", "FontSize"), FontSize=font_size)
set(findall(fig, Type="text"), Color="k")
axes_handles = findall(fig, Type="axes");
for ax = axes_handles'
    ax.XColor = "k";
    ax.YColor = "k";
    ax.ZColor = "k";
    ax.GridColor = [0.7 0.7 0.7];
    ax.GridAlpha = 0.18;
    ax.LabelFontSizeMultiplier = 1;
    ax.Title.FontSize = font_size + 1;
    ax.Title.FontWeight = "normal";
end
set(findall(fig, Type="colorbar"), Color="k")
end
