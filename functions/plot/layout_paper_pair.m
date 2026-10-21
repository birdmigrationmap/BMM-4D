function layout_paper_pair(fig)
% Reserve fixed margins for two appendix panels and their optional colorbars.
fig.Position(4)=1100/300;
axes_handles=findall(fig,Type="axes");
[~,order]=sort(arrayfun(@(ax)ax.Position(1),axes_handles));
axes_handles=axes_handles(order);
for i_panel=1:2
    axes_handles(i_panel).Parent=fig;
    axes_handles(i_panel).Units="normalized";
    axes_handles(i_panel).Position=[0.09+(i_panel-1)*0.50 0.18 0.31 0.72];
    axes_handles(i_panel).XTickLabelRotation=0;
end
bars=findall(fig,Type="colorbar");
[~,order]=sort(arrayfun(@(cb)cb.Position(1),bars));bars=bars(order);
for i_bar=1:numel(bars)
    bars(i_bar).Units="normalized";
    bars(i_bar).Position=[0.415+(i_bar-1)*0.50 0.18 0.012 0.72];
end
end
