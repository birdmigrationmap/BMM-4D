function radar_selected = select_radars_in_bbox(radar_csv, bbox_lonlat)
%SELECT_RADARS_IN_BBOX Read radar metadata and keep sites inside lon/lat bounds.

radar = readtable(radar_csv, TextType="string");
radar_selected = radar( ...
    radar.lon >= bbox_lonlat(1) & radar.lon <= bbox_lonlat(2) & ...
    radar.lat >= bbox_lonlat(3) & radar.lat <= bbox_lonlat(4), :);
end
