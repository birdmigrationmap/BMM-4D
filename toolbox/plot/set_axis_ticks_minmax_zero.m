function set_axis_ticks_minmax_zero(ax, axis_name, limits)
%SET_AXIS_TICKS_MINMAX_ZERO Show only min, max, and 0 ticks when in range.

ticks = unique([limits(1) 0 limits(2)]);
ticks = ticks(ticks >= limits(1) & ticks <= limits(2));

switch lower(axis_name)
    case "x"
        xticks(ax, ticks)
    case "y"
        yticks(ax, ticks)
    case "z"
        zticks(ax, ticks)
    otherwise
        error("set_axis_ticks_minmax_zero:UnknownAxis", ...
            "axis_name must be 'x', 'y', or 'z'.");
end
end
