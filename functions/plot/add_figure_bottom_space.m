function add_figure_bottom_space(fig, extra_px)
% Add room for shared legends without shrinking the existing data panels.
old_height = fig.Position(4) * 300;
fig.Position(4) = (old_height + extra_px) / 300;
objects = [findall(fig, Type="axes"); findall(fig, Type="textboxshape"); findall(fig, Type="rectangleshape")];
for i_object = 1:numel(objects)
    objects(i_object).Units = "normalized";
    position = objects(i_object).Position;
    position(2) = (position(2) * old_height + extra_px) / (old_height + extra_px);
    position(4) = position(4) * old_height / (old_height + extra_px);
    objects(i_object).Position = position;
end
end
